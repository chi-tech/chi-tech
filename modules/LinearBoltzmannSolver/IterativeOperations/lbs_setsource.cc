#include "../lbs_linear_boltzmann_solver.h"

#include "LinearBoltzmannSolver/Groupset/lbs_groupset.h"

#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

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

  const bool apply_fixed_src       = (source_flags & APPLY_FIXED_SOURCES);
  const bool apply_wgs_scatter_src = (source_flags & APPLY_WGS_SCATTER_SOURCES);
  const bool apply_ags_scatter_src = (source_flags & APPLY_AGS_SCATTER_SOURCES);
  const bool apply_wgs_fission_src = (source_flags & APPLY_WGS_FISSION_SOURCES);
  const bool apply_ags_fission_src = (source_flags & APPLY_AGS_FISSION_SOURCES);

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
    if (P0_src && apply_fixed_src)
      src = P0_src->source_value_g.data();

    //======================================== Loop over nodes
    const int num_nodes = transport_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      //=================================== Loop over moments
      for (int m = 0; m < static_cast<int>(num_moments); ++m)
      {
        unsigned int ell = m_to_ell_em_map[m].ell;

        size_t uk_map = transport_view.MapDOF(i, m, 0); //unknown map

        //============================= Loop over groupset groups
        for (size_t g = gs_i; g <= gs_f; ++g)
        {
          double rhs = 0.0;

          //============================== Apply fixed sources
          if (!options.use_src_moments) //using regular material src
            rhs += apply_fixed_src && ell == 0? src[g] : 0.0;
          else if (apply_fixed_src)  //using ext_src_moments
            rhs += ext_src_moments_local[uk_map + g];

          //============================== Apply scattering sources
          const bool moment_avail = ell < S.size();

          //==================== Across groupset
          if (moment_avail && apply_ags_scatter_src)
            for (const auto& [_, gp, sigma_sm] : S[ell].Row(g))
              if (gp < gs_i || gp > gs_f)
                rhs += sigma_sm * phi_old_local[uk_map + gp];

          //==================== Within groupset
          if (moment_avail && apply_wgs_scatter_src)
            for (const auto& [_, gp, sigma_sm] : S[ell].Row(g))
              if (gp >= gs_i && gp <= gs_f)
                rhs += sigma_sm * phi_old_local[uk_map + gp];

          //============================== Apply fission sources
          const bool fission_avail = xs.is_fissionable && ell == 0;

          //==================== Across groupset
          if (fission_avail && apply_ags_fission_src)
          {
            const auto& prod = xs.production_matrix[g];
            for (size_t gp = first_grp; gp <= last_grp; ++gp)
              if (gp < gs_i || gp > gs_f)
              {
                rhs += prod[gp] * phi_old_local[uk_map + gp];

                if (options.use_precursors)
                  for (const auto& precursor : xs.precursors)
                    rhs += precursor.emission_spectrum[g] *
                           precursor.fractional_yield *
                           xs.nu_delayed_sigma_f[gp] *
                           phi_old_local[uk_map + gp];
              }
          }

          //==================== Within groupset
          if (fission_avail && apply_wgs_fission_src)
          {
            const auto& prod = xs.production_matrix[g];
            for (size_t gp = gs_i; gp <= gs_f; ++gp)
            {
              rhs += prod[gp] * phi_old_local[uk_map + gp];

              if (options.use_precursors)
                for (const auto& precursor : xs.precursors)
                  rhs += precursor.emission_spectrum[g] *
                         precursor.fractional_yield *
                         xs.nu_delayed_sigma_f[gp] *
                         phi_old_local[uk_map + gp];
            }
          }

          //============================== Add to destination vector
          destination_q[uk_map + g] += rhs;

        }//for g
      }//for m
    }//for dof i
  }//for cell

  //================================================== Apply point sources
  if (!options.use_src_moments && apply_fixed_src)
    for (const auto& point_source : point_sources)
    {
      const auto& info_list = point_source.ContainingCellsInfo();
      for (const auto& info : info_list)
      {
        auto& transport_view = cell_transport_views[info.cell_local_id];

        const auto& strength = point_source.Strength();
        const auto& node_weights = info.node_weights;
        const double vol_w = info.volume_weight;

        const int num_nodes = transport_view.NumNodes();
        for (int i = 0; i < num_nodes; ++i)
        {
          const size_t uk_map = transport_view.MapDOF(i, /*moment=*/0, /*grp=*/0);
          for (size_t g = gs_i; g <= gs_f; ++g)
            destination_q[uk_map + g] += strength[g] * node_weights[i] * vol_w;
        }//for node i
      }//for cell
    }//for point source

  chi::log.LogEvent(source_event_tag, chi_objects::ChiLog::EventType::EVENT_END);
}
