#include "../k_eigenvalue_solver.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "ChiTimer/chi_timer.h"

#include "chi_mpi.h"
#include "chi_log.h"

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

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
void KEigenvalue::Solver::
SetKSource(LBSGroupset& groupset,
           std::vector<double>& destination_q,
           SourceFlags source_flags)
{
  chi_log.LogEvent(source_event_tag, ChiLog::EventType::EVENT_BEGIN);

  const bool apply_wgs_scatter_src = (source_flags & APPLY_WGS_SCATTER_SOURCE);
  const bool apply_ags_scatter_src = (source_flags & APPLY_AGS_SCATTER_SOURCE);
  const bool apply_wgs_fission_src = (source_flags & APPLY_WGS_FISSION_SOURCE);
  const bool apply_ags_fission_src = (source_flags & APPLY_AGS_FISSION_SOURCE);

  //============================== Get group setup
  int gs_i = groupset.groups[0].id;
  int gs_f = groupset.groups.back().id;

  int first_grp = groups.front().id;
  int last_grp = groups.back().id;

  const auto& m_to_ell_em_map = groupset.quadrature->GetMomentToHarmonicsIndexMap();

  std::vector<double> default_zero_src(groups.size(), 0.0);

  //============================== Loop over local cells
  for (auto& cell : grid->local_cells)
  {
    auto& full_cell_view = cell_transport_views[cell.local_id];

    //==================== Obtain cross-section and src
    int cell_matid = cell.material_id;
    int xs_id = matid_to_xs_map[cell_matid];

    if ((xs_id < 0) || (xs_id >= material_xs.size()))
    {
      chi_log.Log(LOG_ALLERROR)
          << "Cross-section lookup error\n";
      exit(EXIT_FAILURE);
    }

    auto xs = material_xs[xs_id];

    //============================== Loop over nodes
    int num_nodes = full_cell_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      //============================== Loop over moments
      for (int m = 0; m < num_moments; ++m)
      {
        unsigned int ell = m_to_ell_em_map[m].ell;

        size_t ir = full_cell_view.MapDOF(i, m, 0);

        //============================== Loop over groupset groups
        for (size_t g = gs_i; g <= gs_f; ++g)
        {
          //======================================== Apply scattering
          double inscatter_g = 0.0;
          if (ell < xs->transfer_matrices.size())
          {
            //============================== Across-groupset
            if (apply_ags_scatter_src)
            {
              size_t num_transfers =
                  xs->transfer_matrices[ell].rowI_indices[g].size();

              //============================== Loop over transfers
              for (size_t t = 0; t < num_transfers; ++t)
              {
                size_t gprime = xs->transfer_matrices[ell].rowI_indices[g][t];

                if ((gprime < gs_i) or (gprime > gs_f))
                {
                  double sigma_sm = xs->transfer_matrices[ell].rowI_values[g][t];
                  inscatter_g += sigma_sm * phi_old_local[ir + gprime];
                }
              }
            }

            //============================== Within-groupset
            if (apply_wgs_scatter_src)
            {
              size_t num_transfers =
                  xs->transfer_matrices[ell].rowI_indices[g].size();

              //============================== Loop over transfers
              for (size_t t = 0; t < num_transfers; ++t)
              {
                size_t gprime = xs->transfer_matrices[ell].rowI_indices[g][t];

                if ((gprime >= gs_i) and (gprime <= gs_f))
                {
                  double sigma_sm = xs->transfer_matrices[ell].rowI_values[g][t];
                  inscatter_g += sigma_sm * phi_old_local[ir + gprime];
                }
              }
            }
          }//if moment avail
          destination_q[ir + g] += inscatter_g;

          //======================================== Apply fission
          if (xs->is_fissile and (ell == 0))
          {
            double fission_g = 0.0;
            //============================== Across-groupset
            if (apply_ags_fission_src)
            {
              //============================== Loop over groups
              for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
              {
                double nu_sig_f = (options.use_precursors) ?
                                  xs->nu_prompt_sigma_f[gprime] :
                                  xs->nu_sigma_f[gprime];

                if ((gprime < gs_i) or (gprime > gs_f))
                {
                  fission_g += xs->chi[g] * nu_sig_f *
                               phi_prev_local[ir + gprime] / k_eff;

                  //============================== Delayed contributions
                  if (options.use_precursors and xs->num_precursors > 0)
                  {
                    //============================== Loop over precursors
                    for (size_t j = 0; j < xs->num_precursors; ++j)
                    {
                      fission_g += xs->chi_delayed[g][j] *
                                   xs->precursor_yield[j] *
                                   xs->nu_delayed_sigma_f[gprime] *
                                   phi_prev_local[ir + gprime] / k_eff;
                    }
                  }//if use precursors and has precursors
                }//if across groupset
              }//for gprime
            }//if across-groupset

            //============================== Across-groupset
            if (apply_wgs_fission_src)
            {
              //============================== Loop over groups
              for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
              {
                double nu_sig_f = (options.use_precursors) ?
                                  xs->nu_prompt_sigma_f[gprime] : xs->nu_sigma_f[gprime];

                if ((gprime >= gs_i) and (gprime <= gs_f))
                {
                  fission_g += xs->chi[g] * nu_sig_f *
                               phi_prev_local[ir + gprime] / k_eff;

                  //============================== Delayed contributions
                  if (options.use_precursors and xs->num_precursors > 0)
                  {
                    //============================== Loop over precursors
                    for (size_t j = 0; j < xs->num_precursors; ++j)
                    {
                      fission_g += xs->chi_delayed[g][j] *
                                   xs->precursor_yield[j] *
                                   xs->nu_delayed_sigma_f[gprime] *
                                   phi_prev_local[ir + gprime] / k_eff;
                    }
                  }//if use precursors and has precursors
                }//if across groupset
              }//for gprime
            }//if within-groupset
            destination_q[ir + g] += fission_g;

          }//if fissile and ell == 0

        }//for g
      }//for m
    }//for i
  }//for cell

  chi_log.LogEvent(source_event_tag, ChiLog::EventType::EVENT_END);
}