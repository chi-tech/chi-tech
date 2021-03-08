#include "../k_eigenvalue_solver.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "ChiTimer/chi_timer.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

extern double chi_global_timings[20];

using namespace LinearBoltzmann;

//###################################################################
/**Sets the source moments for the groups in the current group set.
 *
 * \param groupset_num Identifies the groupset under consideration.
 * \param apply_mat_src Flag indicating whether the material source needs
 *        to be applied to the routine. This is useful for GMRES
 *        since the material source only features during the computing b.
 *        On this note we also need to treat inscattering this way.
 * \param suppress_phi_old Flag indicating whether to suppress phi_old.
 *
 * */
void KEigenvalue::Solver::SetKSource(LBSGroupset& groupset,
                                     bool apply_mat_src,
                                     bool suppress_phi_old)
{
  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;
  auto pwl = std::static_pointer_cast<SpatialDiscretization_PWLD>(discretization);
  
  bool OneD_Slab = false;

  if (options.geometry_type == GeometryType::ONED_SLAB)
    OneD_Slab = true;

  chi_log.LogEvent(source_event_tag,ChiLog::EventType::EVENT_BEGIN);

  // ----- Groupset group information
  int gs_i = groupset.groups[0].id;
  int gs_f = groupset.groups.back().id;

  int first_grp = groups.front().id;
  int last_grp = groups.back().id;

  std::vector<double> default_zero_src(groups.size(),0.0);

  // ----- Reset source moments
  q_moments_local.assign(q_moments_local.size(),0.0);


  // ----- Loop over local cells
  for (auto& cell : grid->local_cells)
  {
    auto& full_cell_view = cell_transport_views[cell.local_id];

    // ----- Obtain cross-section and src
    int cell_matid = cell.material_id;
    int xs_id = matid_to_xs_map[cell_matid];

    if ((xs_id<0) || (xs_id>=material_xs.size()))
    {
      chi_log.Log(LOG_ALLERROR)
      << "Cross-section lookup error\n";
      exit(EXIT_FAILURE);
    }

    auto xs = material_xs[xs_id];
    
    double inscat_g = 0.0;
    double sigma_sm = 0.0;
    double fission_g = 0.0;
    double precursor_g = 0.0;
    double* q_mom;
    double* phi_oldp;
    double* phi_prevp;
    int gprime;

    // ----- Loop over dofs
    int num_dofs = full_cell_view.dofs;
    for (int i=0; i<num_dofs; i++)
    {
      // ----- Loop over moments
      int m=-1;
      for (int ell=0; ell<=options.scattering_order; ell++)
      {
        int ellmin = -ell;
        int ellmax =  ell;
        if (OneD_Slab)
          {ellmin = 0;ellmax = 0;}

        // ----- Loop over em
        for (int em=ellmin; em<=ellmax; em++)
        {
          m++;
          int ir    = full_cell_view.MapDOF(i,m,0);
          q_mom     = &q_moments_local[ir];
          phi_oldp  = &phi_old_local[ir];
          phi_prevp = &phi_prev_local[ir];
          
          // ----- Loop over groupset groups
          for (int g=gs_i; g<=gs_f; g++)
          {
            // ----- Contribute scattering
            inscat_g = 0.0;
            if ((ell < xs->transfer_matrix.size()) && (!suppress_phi_old) )
            {
              int num_transfers = xs->transfer_matrix[ell].rowI_indices[g].size();
              for (int t=0; t<num_transfers; t++)
              {
                gprime    = xs->transfer_matrix[ell].rowI_indices[g][t];
                sigma_sm  = xs->transfer_matrix[ell].rowI_values[g][t];
                inscat_g += sigma_sm * phi_oldp[gprime];
              }
            }
            q_mom[g] += inscat_g;

            // ----- Contribute fission
            fission_g = 0.0;
            if ((ell == 0) and (apply_mat_src))
            {
              if (xs->is_fissile)
              {                
                for (gprime=first_grp; gprime<=last_grp; ++gprime)
                  if (options.use_precursors)
                    fission_g += xs->chi_g[g]*xs->nu_p_sigma_fg[gprime]*
                                phi_prevp[gprime]/k_eff;
                  else
                    fission_g += xs->chi_g[g]*xs->nu_sigma_fg[gprime]*
                                 phi_prevp[gprime]/k_eff;
              }
            }
            q_mom[g] += fission_g;

            // ----- Contribute precursors
            precursor_g = 0.0;
            if ((ell == 0) and (options.use_precursors))
            {
              if ((apply_mat_src) and (xs->J > 0))
              {
                for (int j=0; j<num_precursors; ++j)
                {
                  for (gprime=first_grp; gprime<=last_grp; ++gprime)
                  {
                    precursor_g += xs->chi_d[g][j]*xs->gamma[j]* 
                                   xs->nu_d_sigma_fg[gprime]*
                                   phi_prevp[gprime]/k_eff;
                  }
                }
              }
            }
            q_mom[g] += precursor_g;
            
          }//for g
        }//for em
      }//for moment
    }//for i
  }//for cell

  chi_log.LogEvent(source_event_tag,ChiLog::EventType::EVENT_END);
}