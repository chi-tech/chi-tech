#include "CHI_MODULES/CHI_NPTRANSPORT/lbs_linear_boltzman_solver.h"

#include "CHI_MODULES/CHI_NPTRANSPORT/SweepChunks/lbs_sweepchunk_pwl_polyhedron.h"
#include "CHI_MESH/CHI_SWEEP/chi_SPDS.h"
#include "CHI_MATH/CHI_DISCRETIZATION/CHI_DISCRETIZATION_PWL/pwl.h"

#include <CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h>
#include <CHI_MESH/CHI_VOLUMEMESHER/Linemesh1D/volmesher_linemesh1d.h>
#include <CHI_MESH/CHI_VOLUMEMESHER/Extruder/volmesher_extruder.h>
#include <CHI_MESH/CHI_VOLUMEMESHER/Predefined2D/volmesher_predefined2d.h>

#include "CHI_TIMER/chi_timer.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern CHI_MPI chi_mpi;
extern CHI_LOG chi_log;

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
void LinearBoltzmanSolver::SetSource(int group_set_num,
                                bool apply_mat_src,
                                bool suppress_phi_old)
{
  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  bool OneD_Slab = false;
  bool TwoD      = false;

  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
    OneD_Slab = true;


  CHI_TIMER t18_setsrctime; t18_setsrctime.Reset();

  //================================================== Get reference to groupset
  NPT_GROUPSET* groupset = group_sets[group_set_num];

  int gs_i = groupset->groups[0]->id;
  int gs_f = groupset->groups.back()->id;

  std::vector<double> default_zero_src(groups.size(),0.0);

  //================================================== Reset source moments
  q_moments_local.assign(q_moments_local.size(),0.0);


  //================================================== Loop over local cells
  for (int c=0; c<grid->local_cell_glob_indices.size(); c++)
  {
    int cell_g_index = grid->local_cell_glob_indices[c];
    auto cell = grid->cells[cell_g_index];

    NPT_CELLVIEW_FULL* full_cell_view =
      (NPT_CELLVIEW_FULL*)cell_transport_views[cell->cell_local_id];

    //=========================================== Obtain cross-section and src
    int cell_matid = cell->material_id;

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
    {
      src = material_srcs[src_id]->source_value_g.data();
    }


    //=========================================== Loop over dofs
    double inscat_g = 0.0;
    double sigma_sm = 0.0;
    double* q_mom;
    double* phi_old;
    int num_dofs = full_cell_view->dofs;
    int gprime;
    for (int i=0; i<num_dofs; i++)
    {
      //==================================== Loop over moments
      int m=-1;
      for (int ell=0; ell<=options.scattering_order; ell++)
      {
        int ellmin = -ell;
        int ellmax =  ell;
        if (OneD_Slab)
          {ellmin = 0;ellmax = 0;}

        for (int em=ellmin; em<=ellmax; em++)
        {
          m++;
          int ir = full_cell_view->MapDOF(i,m,0);
          q_mom   = &q_moments_local.data()[ir];
          phi_old = &phi_old_local.data()[ir];

          //============================= Loop over groupset groups
          for (int g=gs_i; g<=gs_f; g++)
          {
            if (apply_mat_src && (m==0))
              q_mom[g] += src[g];


            inscat_g = 0.0;
            //====================== Apply across-groupset scattering
            if ((ell < xs->transfer_matrix.size()) && (apply_mat_src) )
            {
              int num_transfers = xs->transfer_matrix[ell].inds_rowI[g].size();
              for (int t=0; t<num_transfers; t++)
              {
                gprime    = xs->transfer_matrix[ell].inds_rowI[g][t];
                if ((gprime < gs_i) || (gprime > gs_f))
                {
                  sigma_sm  = xs->transfer_matrix[ell].rowI_colJ[g][t];
                  inscat_g += sigma_sm*phi_old[gprime];
                }
              }
            }//if moment avail

            //====================== Apply within-groupset scattering
            if ((ell < xs->transfer_matrix.size()) && (!suppress_phi_old) )
            {
              int num_transfers = xs->transfer_matrix[ell].inds_rowI[g].size();
              for (int t=0; t<num_transfers; t++)
              {
                gprime    = xs->transfer_matrix[ell].inds_rowI[g][t];
                if ((gprime >= gs_i) && (gprime<=gs_f))
                {
                  sigma_sm  = xs->transfer_matrix[ell].rowI_colJ[g][t];
                  inscat_g += sigma_sm*phi_old[gprime];
                }
              }
            }//if moment avail
            q_mom[g] += inscat_g;
          }//for g
        }

      }//for moment
    }//for dof i

  }//for cell

  chi_global_timings[18] += t18_setsrctime.GetTime();
  chi_global_timings[19] += 1.0;
}