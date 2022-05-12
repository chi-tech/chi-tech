#include "../lbs_linear_boltzmann_solver.h"

#include "ChiTimer/chi_timer.h"

#include <chi_mpi.h>
#include <chi_log.h>


;

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
void lbs::SteadySolver::
  SetSource(LBSGroupset& groupset,
            std::vector<double>& destination_q,
            SourceFlags source_flags)
{
  chi::log.LogEvent(source_event_tag, chi_objects::ChiLog::EventType::EVENT_BEGIN);

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
  // Apply all nodal sources
  for (const auto& cell : grid->local_cells)
  {
    auto& transport_view = cell_transport_views[cell.local_id];

    //==================== Obtain xs
    auto xs = transport_view.XS();
    auto P0_src = matid_to_src_map[cell.material_id];

    const auto& S = xs.transfer_matrices;

    //==================== Obtain src
    double* src = default_zero_src.data();
    if (P0_src and apply_mat_src)
      src = P0_src->source_value_g.data();

    //=========================================== Loop over nodes
    const int num_nodes = transport_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      //==================================== Loop over moments
      for (int m = 0; m < static_cast<int>(num_moments); ++m)
      {
        unsigned int ell = m_to_ell_em_map[m].ell;

        size_t uk_map = transport_view.MapDOF(i, m, 0); //unknown map

        //============================= Loop over groupset groups
        for (size_t g = gs_i; g <= gs_f; ++g)
        {
          if (not options.use_src_moments) //using regular material src
          {
            if (apply_mat_src and ell == 0)
              destination_q[uk_map + g] += src[g];
          }
          else if (apply_mat_src)  //using ext_src_moments
            destination_q[uk_map + g] += ext_src_moments_local[uk_map + g];

          double inscatter_g = 0.0;
          const bool moment_avail = (ell < S.size());

          //====================== Apply across-groupset scattering
          if (moment_avail and apply_ags_scatter_src)
          {
            size_t num_transfers =
                S[ell].rowI_indices[g].size();

            //=============== Loop over transfers
            for (size_t t = 0; t < num_transfers; ++t)
            {
              size_t gprime =
                  S[ell].rowI_indices[g][t];

              if ((gprime < gs_i) or (gprime > gs_f))
              {
                double sigma_sm = S[ell].rowI_values[g][t];
                inscatter_g += sigma_sm * phi_old_local[uk_map + gprime];
              }
            }
          }//if moment_avail

          //====================== Apply within-groupset scattering
          if (moment_avail and apply_wgs_scatter_src)
          {
            size_t num_transfers =
                S[ell].rowI_indices[g].size();

            //=============== Loop over transfers
            for (size_t t = 0; t < num_transfers; ++t)
            {
              size_t gprime = S[ell].rowI_indices[g][t];
              if ((gprime >= gs_i) and (gprime <= gs_f))
              {
                double sigma_sm = S[ell].rowI_values[g][t];
                inscatter_g += sigma_sm * phi_old_local[uk_map + gprime];
              }
            }
          }
          destination_q[uk_map + g] += inscatter_g;


          double infission_g = 0.0;
          const bool fission_avail = (xs.is_fissile and ell == 0);

          //====================== Apply accross-groupset fission
          if (fission_avail and apply_ags_fission_src)
          {
            //=============== Loop over groups
            for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
            {
              if ((gprime < gs_i) or (gprime > gs_f))
              {
                //without delayed neutron precursors
                if (not options.use_precursors)
                  infission_g += xs.chi[g] *
                                 xs.nu_sigma_f[gprime] *
                                 phi_old_local[uk_map + gprime];

                //with delayed neutron precursors
                else
                {
                  //Prompt fission
                  infission_g += xs.chi_prompt[g] *
                                 xs.nu_prompt_sigma_f[gprime] *
                                 phi_old_local[uk_map + gprime];

                  //Delayed fission
                  for (size_t j = 0; j < xs.num_precursors; ++j)
                    infission_g += xs.chi_delayed[g][j] *
                                   xs.precursor_yield[j] *
                                   xs.nu_delayed_sigma_f[gprime] *
                                   phi_old_local[uk_map + gprime];
                }
              }
            }//for gprime
          }//if zeroth moment

          //====================== Apply within-groupset fission
          if (fission_avail and apply_wgs_fission_src)
          {
            //=============== Loop over groups
            for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
            {
              if ((gprime >= gs_i) and (gprime <= gs_f))
              {
                //without delayed neutron precursors
                if (not options.use_precursors)
                  infission_g += xs.chi[g] *
                                 xs.nu_sigma_f[gprime] *
                                 phi_old_local[uk_map + gprime];

                //with delayed neutron precursors
                else
                {
                  //Prompt fission
                  infission_g += xs.chi_prompt[g] *
                                 xs.nu_prompt_sigma_f[gprime] *
                                 phi_old_local[uk_map + gprime];

                  //Delayed fission
                  for (size_t j = 0; j < xs.num_precursors; ++j)
                    infission_g += xs.chi_delayed[g][j] *
                                   xs.precursor_yield[j] *
                                   xs.nu_delayed_sigma_f[gprime] *
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

  //================================================== Apply point sources
  if ((not options.use_src_moments) and apply_mat_src)
  {
    for (const auto& point_source : point_sources)
    {
      if (not point_source.LocallyOwned()) continue;
      const uint64_t cell_local_id = point_source.OwningCellLocalID();

      auto& transport_view = cell_transport_views[cell_local_id];

      const auto& strength = point_source.Strength();
      const auto& node_weights = point_source.NodeWeights();

      const int num_nodes = transport_view.NumNodes();
      for (int i = 0; i < num_nodes; ++i)
      {
        const size_t uk_map = transport_view.MapDOF(i, /*moment=*/0, /*grp=*/0);
        for (size_t g = gs_i; g <= gs_f; ++g)
          destination_q[uk_map + g] += strength[g] * node_weights[i];
      }//for node i
    }//for point source
  }//if apply mat src

  chi::log.LogEvent(source_event_tag, chi_objects::ChiLog::EventType::EVENT_END);
}
