#include "../lbts_transient_solver.h"

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
void lbs::TransientSolver::
  SetTransientSource(LBSGroupset& groupset,
                     std::vector<double>& destination_q,
                     SourceFlags source_flags)
{
  chi::log.LogEvent(source_event_tag, chi_objects::ChiLog::EventType::EVENT_BEGIN);

  const auto& BackwardEuler = chi_math::SteppingMethod::BACKWARD_EULER;
  const auto& CrankNicolson = chi_math::SteppingMethod::CRANK_NICHOLSON;

  double theta;
  if      (method == BackwardEuler) theta = 1.0;
  else if (method == CrankNicolson) theta = 0.5;
  else                              theta = 0.7;

  const double eff_dt = theta * dt;

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
    const auto& fe_values = unit_cell_matrices[cell.local_id];
    auto& transport_view = cell_transport_views[cell.local_id];
    const double cell_volume = transport_view.Volume();

    //==================== Obtain xs
    auto xs = transport_view.XS();
    auto P0_src = matid_to_src_map[cell.material_id];

    const auto& S = xs.transfer_matrices;

    //==================== Obtain src
    double* src = default_zero_src.data();
    if (P0_src and apply_fixed_src)
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
            if (apply_fixed_src and ell == 0)
              destination_q[uk_map + g] += src[g];
          }
          else if (apply_fixed_src)  //using ext_src_moments
            destination_q[uk_map + g] += ext_src_moments_local[uk_map + g];

          double inscatter_g = 0.0;
          const bool moment_avail = (ell < S.size());

          //====================== Apply across-groupset scattering
          if (moment_avail and apply_ags_scatter_src)
            for (const auto& [row_g, gprime, sigma_sm] : S[ell].Row(g))
              if ((gprime < gs_i) or (gprime > gs_f))
                inscatter_g += sigma_sm * phi_old_local[uk_map + gprime];

          //====================== Apply within-groupset scattering
          if (moment_avail and apply_wgs_scatter_src)
            for (const auto& [row_g, gprime, sigma_sm] : S[ell].Row(g))
              if ((gprime >= gs_i) and (gprime <= gs_f))
                inscatter_g += sigma_sm * phi_old_local[uk_map + gprime];

          destination_q[uk_map + g] += inscatter_g;

          double infission_g = 0.0;
          const bool fission_avail = (xs.is_fissionable and ell == 0);

          //====================== Apply accross-groupset fission
          if (fission_avail and apply_ags_fission_src)
          {
            const double spec = options.use_precursors?
                xs.chi_prompt[g] : xs.chi[g];

            //=============== Loop over groups
            for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
              if ((gprime < gs_i) or (gprime > gs_f))
              {
                const double nu_sigma_f = (options.use_precursors)?
                  xs.nu_prompt_sigma_f[gprime] : xs.nu_sigma_f[gprime];

                  infission_g += spec * nu_sigma_f *
                                 phi_old_local[uk_map + gprime];
              }//if gprime outside current groupset

            //=================================== Delayed fission
            if (options.use_precursors)
            {
              for (const auto& precursor : xs.precursors)
              {
                double coeff =
                    precursor.emission_spectrum[g] *
                    precursor.decay_constant /
                    (1.0 + eff_dt * precursor.decay_constant);

                for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
                {
                  if ((gprime < gs_i) or (gprime > gs_f))
                    infission_g += coeff * eff_dt *
                                   precursor.fractional_yield *
                                   xs.nu_delayed_sigma_f[gprime] *
                                   phi_old_local[uk_map + gprime] /
                                   cell_volume;
                }//for gprime
              }//for precursors
            }//if use precursors
          }//if fission_avail and apply_ags_fission_src

          //====================== Apply within-groupset fission
          if (fission_avail and apply_wgs_fission_src)
          {
            const double spec = options.use_precursors?
                xs.chi_prompt[g] : xs.chi[g];

            //=============== Loop over groups apply promp fission
            for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
              if ((gprime >= gs_i) and (gprime <= gs_f))
              {
                const double nu_sigma_f = (options.use_precursors)?
                  xs.nu_prompt_sigma_f[gprime] : xs.nu_sigma_f[gprime];

                infission_g += spec * nu_sigma_f *
                               phi_old_local[uk_map + gprime];
              }//if gprime inside current groupset

            //=================================== Apply Delayed fission
            if (options.use_precursors)
            {
              for (const auto& precursor : xs.precursors)
              {
                const double coeff =
                    precursor.emission_spectrum[g] *
                    precursor.decay_constant /
                    (1.0 + eff_dt * precursor.decay_constant);

                //==================== Delayed fission rate contributions
                for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
                {
                  if ((gprime >= gs_i) and (gprime <= gs_f))
                    infission_g += coeff * eff_dt *
                                   precursor.fractional_yield *
                                   xs.nu_delayed_sigma_f[gprime] *
                                   phi_old_local[uk_map + gprime] /
                                   cell_volume;
                }//for gprime
              }//for precursors
            }//if use precursors
          }//if fission_avail and apply_wgs_fission_src

          //====================== Apply precursors
          if (fission_avail and apply_fixed_src and options.use_precursors)
          {
            const auto& J = max_precursors_per_material;
            for (unsigned int j = 0; j < xs.num_precursors; ++j)
            {
              const size_t dof_map = cell.local_id * J + j;
              const auto& precursor = xs.precursors[j];

              const double coeff =
                  precursor.emission_spectrum[g] *
                  precursor.decay_constant /
                  (1.0 + eff_dt * precursor.decay_constant);

              infission_g += coeff * precursor_prev_local[dof_map];
            }//for precursors
          }//if fission_avail and apply_wgs_fission_src

          destination_q[uk_map + g] += infission_g;

        }//for g
      }//for m
    }//for dof i
  }//for cell

  //================================================== Apply point sources
  if ((not options.use_src_moments) and apply_fixed_src)
  {
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
  }//if apply mat src

  chi::log.LogEvent(source_event_tag, chi_objects::ChiLog::EventType::EVENT_END);
}