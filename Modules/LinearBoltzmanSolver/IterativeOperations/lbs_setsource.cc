#include "../lbs_linear_boltzman_solver.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>

#include "ChiTimer/chi_timer.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

extern double chi_global_timings[20];

//###################################################################
/**Sets the source moments for the groups in the current group set.
 *
 * \param group_set_num Identifies the groupset under consideration.
 * \param apply_mat_src Flag indicating whether the material source needs
 *        to be applied to the routine. This is useful for GMRES
 *        since the material source only features during the computing b.
 *        On this note we also need to treat inscattering this way.
 * \param suppress_phi_old Flag indicating whether to suppress phi_old.
 *
 * */
void LinearBoltzman::Solver::SetSource(int group_set_num,
                                bool apply_mat_src,
                                bool suppress_phi_old)
{
  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  bool OneD_Slab = false;

  if (options.geometry_type == GeometryType::ONED_SLAB)
    OneD_Slab = true;

  chi_log.LogEvent(source_event_tag,ChiLog::EventType::EVENT_BEGIN);

  //================================================== Get reference to groupset
  LBSGroupset* groupset = group_sets[group_set_num];

  int gs_i = groupset->groups[0]->id;
  int gs_f = groupset->groups.back()->id;

  int first_grp = groups.front()->id;
  int last_grp = groups.back()->id;

  std::vector<double> default_zero_src(groups.size(),0.0);

  //================================================== Reset source moments
  q_moments_local.assign(q_moments_local.size(),0.0);


  //================================================== Loop over local cells
  for (const auto& cell : grid->local_cells)
  {
    auto full_cell_view =
      (LinearBoltzman::CellViewFull*)cell_transport_views[cell.local_id];

    //=========================================== Obtain cross-section and src
    int cell_matid = cell.material_id;

    int xs_id = matid_to_xs_map[cell_matid];
    int src_id= matid_to_src_map[cell_matid];

    if ((xs_id<0) || (xs_id>=material_xs.size()))
    {
      chi_log.Log(LOG_ALLERROR)
      << "Cross-section lookup error\n";
      exit(EXIT_FAILURE);
    }

    chi_physics::TransportCrossSections* xs =
      material_xs[xs_id];

    //=========================================== Obtain material source
    double* src = default_zero_src.data();
    if ( (src_id >= 0) && (apply_mat_src) )
      src = material_srcs[src_id]->source_value_g.data();


    //=========================================== Loop over dofs
    double inscat_g = 0.0;
    double sigma_sm = 0.0;
    double* q_mom;
    double* phi_oldp;
    int num_dofs = full_cell_view->dofs;
    int gprime;
    for (int i=0; i<num_dofs; i++)
    {
      //==================================== Loop over moments
      int m=-1;
      for (int ell=0; ell<=options.scattering_order; ell++)
      {
        int ellmin = OneD_Slab? 0 : -ell;
        int ellmax = OneD_Slab? 0 :  ell;

        for (int em=ellmin; em<=ellmax; em++)
        {
          m++;
          int ir = full_cell_view->MapDOF(i,m,0);
          q_mom    = &q_moments_local[ir];
          phi_oldp = &phi_old_local[ir];

          //============================= Loop over groupset groups
          for (int g=gs_i; g<=gs_f; g++)
          {
            if (apply_mat_src && (m==0))
              q_mom[g] += src[g];


            inscat_g = 0.0;
            //====================== Apply across-groupset scattering
            if ((ell < xs->transfer_matrix.size()) && (apply_mat_src) )
            {
              int num_transfers = xs->transfer_matrix[ell].rowI_indices[g].size();
              for (int t=0; t<num_transfers; t++)
              {
                gprime    = xs->transfer_matrix[ell].rowI_indices[g][t];
                if ((gprime < gs_i) || (gprime > gs_f))
                {
                  sigma_sm  = xs->transfer_matrix[ell].rowI_values[g][t];
                  inscat_g += sigma_sm * phi_oldp[gprime];
                }
              }
            }//if moment avail

            //====================== Apply within-groupset scattering
            if ((ell < xs->transfer_matrix.size()) && (!suppress_phi_old) )
            {
              int num_transfers = xs->transfer_matrix[ell].rowI_indices[g].size();
              for (int t=0; t<num_transfers; t++)
              {
                gprime    = xs->transfer_matrix[ell].rowI_indices[g][t];
                if ((gprime >= gs_i) && (gprime<=gs_f))
                {
                  sigma_sm  = xs->transfer_matrix[ell].rowI_values[g][t];
                  inscat_g += sigma_sm * phi_oldp[gprime];
                }
              }
            }//if moment avail

            q_mom[g] += inscat_g;

            //====================== Apply accross-groupset fission
            if ((ell == 0) and (apply_mat_src))
            {
              for (gprime=first_grp; gprime<=last_grp; ++gprime)
              {
                if ((gprime < gs_i) || (gprime > gs_f))
                {
                  q_mom[g] += xs->chi_g[g]*
                              xs->nu_sigma_fg[gprime]*
                              phi_oldp[gprime];
                }
              }//for gprime
            }//if zeroth moment

            //====================== Apply within-groupset fission
            if ((ell == 0) and (!suppress_phi_old))
            {
              for (gprime=first_grp; gprime<=last_grp; ++gprime)
              {
                if ((gprime >= gs_i) && (gprime<=gs_f))
                {
                  q_mom[g] += xs->chi_g[g]*
                              xs->nu_sigma_fg[gprime]*
                              phi_oldp[gprime];
                }
              }//for gprime
            }//if zeroth moment
          }//for g
        }

      }//for moment
    }//for dof i

  }//for cell

  chi_log.LogEvent(source_event_tag,ChiLog::EventType::EVENT_END);
}