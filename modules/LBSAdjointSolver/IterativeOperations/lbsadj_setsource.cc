#include "LBSAdjointSolver/lbsadj_solver.h"

#include "ChiTimer/chi_timer.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "LinearBoltzmannSolver/Groupset/lbs_groupset.h"

void lbs_adjoint::AdjointSolver::
  SetAdjointSource(lbs::LBSGroupset &groupset,
                   std::vector<double> &destination_q,
                   lbs::SourceFlags source_flags)
{
  chi::log.LogEvent(source_event_tag, chi_objects::ChiLog::EventType::EVENT_BEGIN);

  using Flag = lbs::SourceFlags;
  const bool apply_fixed_src       = (source_flags & Flag::APPLY_FIXED_SOURCES);
  const bool apply_wgs_scatter_src = (source_flags & Flag::APPLY_WGS_SCATTER_SOURCES);
  const bool apply_ags_scatter_src = (source_flags & Flag::APPLY_AGS_SCATTER_SOURCES);
  const bool apply_wgs_fission_src = (source_flags & Flag::APPLY_WGS_FISSION_SOURCES);
  const bool apply_ags_fission_src = (source_flags & Flag::APPLY_AGS_FISSION_SOURCES);

  //================================================== Get group setup
  auto gs_i = static_cast<size_t>(groupset.groups[0].id);
  auto gs_f = static_cast<size_t>(groupset.groups.back().id);

  auto first_grp = static_cast<size_t>(groups.front().id);
  auto last_grp = static_cast<size_t>(groups.back().id);

  const auto& m_to_ell_em_map =
    groupset.quadrature->GetMomentToHarmonicsIndexMap();

  std::vector<double> default_zero_src(groups.size(), 0.0);

  //================================================== Loop over local cells
  //                                                   and apply scattering only
  for (const auto& cell : grid->local_cells)
  {
    auto& full_cell_view = cell_transport_views[cell.local_id];

    //==================== Obtain xs
    auto xs     = full_cell_view.XS();
    auto P0_src = matid_to_src_map[cell.material_id];

    const auto& S = matid_to_S_transpose[cell.material_id];

    //=========================================== Loop over nodes
    const int num_nodes = full_cell_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      //==================================== Loop over moments
      for (int m = 0; m < static_cast<int>(num_moments); ++m)
      {
        unsigned int ell = m_to_ell_em_map[m].ell;

        size_t uk_map = full_cell_view.MapDOF(i, m, 0); //unknown map

        //============================= Loop over groupset groups
        for (size_t g = gs_i; g <= gs_f; ++g)
        {
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

  //================================================== Apply reference QOI
  if (apply_fixed_src)
  {
    for (const auto& qoi_data : response_functions)
    {
      const auto& qoi_designation = qoi_data.first;
      const auto& qoi_cell_subscription = qoi_data.second;

      if (qoi_designation.name == basic_options("REFERENCE_RF").StringValue())
      {
        for (size_t local_id : qoi_cell_subscription)
        {
          const auto& full_cell_view = cell_transport_views[local_id];
          const auto& cell = grid->local_cells[local_id];
          const auto& response = qoi_designation.GetMGResponse(cell, num_groups);
          const int num_nodes = full_cell_view.NumNodes();
          for (int i = 0; i < num_nodes; ++i)
          {
            size_t uk_map = full_cell_view.MapDOF(i, 0, 0); //unknown map

            for (size_t g = gs_i; g <= gs_f; ++g)
            {
              destination_q[uk_map + g] += response[g];
            }//for group
          }//for node
        }//for local cell-id of qoi
      }//if ref-qoi
    }//for qoi
  }

  chi::log.LogEvent(source_event_tag, chi_objects::ChiLog::EventType::EVENT_END);
}