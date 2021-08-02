#include "../lbs_linear_boltzmann_solver.h"

#include "ChiTimer/chi_timer.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

//###################################################################
/**Sets the source moments for the groups in the current group set.
 *
 * \param groupset The groupset the under consideration.
 * \param destination_q A vector to contribute the source to.
 * \param source_flags Flags for adding specific terms into the
 *        destination vector. Available flags are for applying
 *        the material source, across/within-group scattering,
 *        and across/within-groups fission.
 *
 * */
void LinearBoltzmann::Solver::
  SetSource(LBSGroupset& groupset,
            std::vector<double>& destination_q,
            SourceFlags source_flags)
{
  chi_log.LogEvent(source_event_tag, ChiLog::EventType::EVENT_BEGIN);

  const bool apply_mat_src         = (source_flags & APPLY_MATERIAL_SOURCE);
  const bool apply_wgs_scatter_src = (source_flags & APPLY_WGS_SCATTER_SOURCE);
  const bool apply_ags_scatter_src = (source_flags & APPLY_AGS_SCATTER_SOURCE);
  const bool apply_wgs_fission_src = (source_flags & APPLY_WGS_FISSION_SOURCE);
  const bool apply_ags_fission_src = (source_flags & APPLY_AGS_FISSION_SOURCE);

  //================================================== Get group setup
  auto gs_i = static_cast<size_t>(groupset.groups[0].id);
  auto gs_f = static_cast<size_t>(groupset.groups.back().id);

  auto first_grp = static_cast<size_t>(groups.front().id);
  auto last_grp = static_cast<size_t>(groups.back().id);

  const auto& m_to_ell_em_map =
      groupset.quadrature->GetMomentToHarmonicsIndexMap();

  std::vector<double> default_zero_src(groups.size(), 0.0);

  //================================================== Loop over local cells
  for (const auto& cell : grid->local_cells)
  {
    auto& full_cell_view = cell_transport_views[cell.local_id];

    //==================== Obtain xs
    int cell_matid = cell.material_id;
    int xs_id = matid_to_xs_map[cell_matid];
    int src_id = matid_to_src_map[cell_matid];

    int num_mat_xs = static_cast<int>(material_xs.size());

    if ((xs_id < 0) || (xs_id >= num_mat_xs))
    {
      chi_log.Log(LOG_ALLERROR)
          << "Cross-section lookup error\n";
      exit(EXIT_FAILURE);
    }

    auto xs = material_xs[xs_id];

    //==================== Obtain src
    double* src = default_zero_src.data();
    if ((src_id >= 0) && (apply_mat_src))
      src = material_srcs[src_id]->source_value_g.data();

    //======================================== Loop over nodes
    const int num_nodes = full_cell_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      //=================================== Loop over moments
      for (int m = 0; m < num_moments; ++m)
      {
        unsigned int ell = m_to_ell_em_map[m].ell;

        size_t uk_map = full_cell_view.MapDOF(i, m, 0); //unknown map

        //============================== Loop over groupset groups
        for (size_t g = gs_i; g <= gs_f; ++g)
        {
          if ( apply_mat_src and (ell == 0) and (not options.use_src_moments))
            destination_q[uk_map + g] += src[g];
          else if (apply_mat_src and options.use_src_moments)
            destination_q[uk_map + g] += ext_src_moments_local[uk_map + g];

          double inscatter_g = 0.0;
          const bool moment_avail = (ell < xs->transfer_matrices.size());

          //=================================== Apply across-groupset scattering
          if (moment_avail and apply_ags_scatter_src)
          {
            size_t num_transfers =
                xs->transfer_matrices[ell].rowI_indices[g].size();

            //============================== Loop over transfers
            for (size_t t = 0; t < num_transfers; ++t)
            {
              size_t gprime =
                  xs->transfer_matrices[ell].rowI_indices[g][t];

              if ((gprime < gs_i) or (gprime > gs_f))
              {
                double sigma_sm = xs->transfer_matrices[ell].rowI_values[g][t];
                inscatter_g += sigma_sm * phi_old_local[uk_map + gprime];
              }
            }
          }//if moment_avail

          //=================================== Apply within-groupset scattering
          if (moment_avail and apply_wgs_scatter_src)
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
                inscatter_g += sigma_sm * phi_old_local[uk_map + gprime];
              }
            }
          }
          destination_q[uk_map + g] += inscatter_g;


          double infission_g = 0.0;
          const bool fission_avail = (xs->is_fissile and ell == 0);

          //=================================== Apply accross-groupset fission
          if (fission_avail and apply_ags_fission_src)
          {
            //================================ Loop over groups
            for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
            {
              if ((gprime < gs_i) or (gprime > gs_f))
              {
                //without delayed neutron precursors
                if (not options.use_precursors)
                  infission_g += xs->chi[g] *
                                 xs->nu_sigma_f[gprime] *
                                 phi_old_local[uk_map + gprime];

                  //with delayed neutron precursors
                else
                {
                  //==================== Prompt fission
                  infission_g += xs->chi_prompt[g] *
                                 xs->nu_prompt_sigma_f[gprime] *
                                 phi_old_local[uk_map + gprime];

                  //==================== Delayed fission
                  for (size_t j = 0; j < xs->num_precursors; ++j)
                    infission_g += xs->chi_delayed[g][j] *
                                   xs->precursor_yield[j] *
                                   xs->nu_delayed_sigma_f[gprime] *
                                   phi_old_local[uk_map + gprime];
                }
              }
            }//for gprime
          }//if zeroth moment

          //=================================== Apply within-groupset fission
          if (fission_avail and apply_wgs_fission_src)
          {
            //============================== Loop over groups
            for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
            {
              if ((gprime >= gs_i) and (gprime <= gs_f))
              {
                //without delayed neutron precursors
                if (not options.use_precursors)
                  infission_g += xs->chi[g] *
                                 xs->nu_sigma_f[gprime] *
                                 phi_old_local[uk_map + gprime];

                  //with delayed neutron precursors
                else
                {
                  //==================== Prompt fission
                  infission_g += xs->chi_prompt[g] *
                                 xs->nu_prompt_sigma_f[gprime] *
                                 phi_old_local[uk_map + gprime];

                  //==================== Delayed fission
                  for (size_t j = 0; j < xs->num_precursors; ++j)
                    infission_g += xs->chi_delayed[g][j] *
                                   xs->precursor_yield[j] *
                                   xs->nu_delayed_sigma_f[gprime] *
                                   phi_old_local[uk_map + gprime];
                }
              }
            }
          }
          destination_q[uk_map + g] += infission_g;

        }//for g
      }//for m
    }//for dof i
  }//for cell

  chi_log.LogEvent(source_event_tag, ChiLog::EventType::EVENT_END);
}