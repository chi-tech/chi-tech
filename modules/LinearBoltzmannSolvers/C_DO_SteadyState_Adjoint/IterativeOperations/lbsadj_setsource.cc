#include "C_DO_SteadyState_Adjoint/lbsadj_solver.h"

#include "ChiTimer/chi_timer.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

void lbs::DiscOrdSteadyStateAdjointSolver::
  SetAdjointSource(lbs::LBSGroupset &groupset,
                   std::vector<double> &destination_q,
                   const std::vector<double>& phi,
                   lbs::SourceFlags source_flags)
{
  chi::log.LogEvent(source_event_tag_, chi_objects::ChiLog::EventType::EVENT_BEGIN);

  using Flag = lbs::SourceFlags;
  const bool apply_fixed_src       = (source_flags & Flag::APPLY_FIXED_SOURCES);
  const bool apply_wgs_scatter_src = (source_flags & Flag::APPLY_WGS_SCATTER_SOURCES);
  const bool apply_ags_scatter_src = (source_flags & Flag::APPLY_AGS_SCATTER_SOURCES);
  const bool apply_wgs_fission_src = (source_flags & Flag::APPLY_WGS_FISSION_SOURCES);
  const bool apply_ags_fission_src = (source_flags & Flag::APPLY_AGS_FISSION_SOURCES);

  //================================================== Get group setup
  auto gs_i = static_cast<size_t>(groupset.groups_.front().id_);
  auto gs_f = static_cast<size_t>(groupset.groups_.back().id_);

  auto first_grp = static_cast<size_t>(groups_.front().id_);
  auto last_grp = static_cast<size_t>(groups_.back().id_);

  const auto& m_to_ell_em_map =
    groupset.quadrature->GetMomentToHarmonicsIndexMap();

  std::vector<double> default_zero_src(groups_.size(), 0.0);

  //================================================== Loop over local cells
  //                                                   and apply scattering only
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& full_cell_view = cell_transport_views_[cell.local_id];

    //==================== Obtain xs
    auto xs     = full_cell_view.XS();
    auto P0_src = matid_to_src_map_[cell.material_id];

    const auto& S = matid_to_S_transpose_[cell.material_id];

    //======================================== Loop over nodes
    const int num_nodes = full_cell_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      //=================================== Loop over moments
      for (int m = 0; m < static_cast<int>(num_moments_); ++m)
      {
        unsigned int ell = m_to_ell_em_map[m].ell;

        size_t uk_map = full_cell_view.MapDOF(i, m, 0); //unknown map

        //============================== Loop over groupset groups
        for (size_t g = gs_i; g <= gs_f; ++g)
        {
          double rhs = 0.0;

          //============================== Apply scattering sources
          const bool moment_avail = ell < S.size();

          //==================== Across groupset
          if (moment_avail and apply_ags_scatter_src)
            for (const auto& [_, gp, sigma_sm] : S[ell].Row(g))
              if (gp < gs_i or gp > gs_f)
                rhs += sigma_sm * phi[uk_map + gp];

          //==================== Within groupset
          if (moment_avail and apply_wgs_scatter_src)
            for (const auto& [_, gp, sigma_sm] : S[ell].Row(g))
              if (gp >= gs_i and gp <= gs_f)
                rhs += sigma_sm * phi[uk_map + gp];

          //============================== Apply fission sources
          const bool fission_avail = xs.is_fissionable and ell == 0;

          //==================== Across groupset
          if (fission_avail and apply_ags_fission_src)
          {
            const auto& prod = xs.production_matrix[g];
            for (size_t gp = first_grp; gp <= last_grp; ++gp)
              if (gp < gs_i or gp > gs_f)
              {
                rhs += prod[gp] * phi[uk_map + gp];

                if (options_.use_precursors)
                  for (const auto& precursor: xs.precursors)
                    rhs += precursor.emission_spectrum[g] *
                           precursor.fractional_yield *
                           xs.nu_delayed_sigma_f[gp] *
                           phi[uk_map + gp];
              }
          }

          //==================== Within groupset
          if (fission_avail and apply_wgs_fission_src)
          {
            const auto& prod = xs.production_matrix[g];
            for (size_t gp = gs_i; gp <= gs_f; ++gp)
            {
              rhs += prod[gp] * phi[uk_map + gp];

              if (options_.use_precursors)
                for (const auto& precursor: xs.precursors)
                  rhs += precursor.emission_spectrum[g] *
                         precursor.fractional_yield *
                         xs.nu_delayed_sigma_f[gp] *
                         phi[uk_map + gp];
            }
          }

          //============================== Add to destination vector
          destination_q[uk_map + g] += rhs;

        }//for g
      }//for m
    }//for dof i
  }//for cell

  //================================================== Apply reference QOI
  if (apply_fixed_src)
    for (const auto& qoi_data : response_functions_)
    {
      const auto& qoi_designation = qoi_data.first;
      const auto& qoi_cell_subscription = qoi_data.second;

      if (qoi_designation.name == basic_options_("REFERENCE_RF").StringValue())
      {
        for (size_t local_id : qoi_cell_subscription)
        {
          const auto& full_cell_view = cell_transport_views_[local_id];
          const auto& cell = grid_ptr_->local_cells[local_id];
          const auto& response = qoi_designation.GetMGResponse(cell, num_groups_);
          const int num_nodes = full_cell_view.NumNodes();
          for (int i = 0; i < num_nodes; ++i)
          {
            size_t uk_map = full_cell_view.MapDOF(i, 0, 0); //unknown map
            for (size_t g = gs_i; g <= gs_f; ++g)
              destination_q[uk_map + g] += response[g];
          }//for node
        }//for local cell-id of qoi
      }//if ref-qoi
    }//for qoi

  chi::log.LogEvent(source_event_tag_, chi_objects::ChiLog::EventType::EVENT_END);
}