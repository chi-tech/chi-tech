#include "../k_eigenvalue_solver.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"


#include "chi_mpi.h"
#include "chi_log.h"

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

using namespace LinearBoltzmann;

//###################################################################
/**Sets the source moments for the groups in the current group set.
 *
 * \param groupset_num Identifies the groupset under consideration.
 * \param destination_q A vector to contribute the source to.
 * \param source_flags Flags for adding specific terms into the
 *        destination vector. Available flags are for applying
 *        the material source, across/within-group scattering,
 *        and across/within-groups fission.
 *
 * */
void KEigenvalue::Solver::SetKSource(LBSGroupset& groupset,
                                     std::vector<double>& destination_q,
                                      SourceFlags source_flags)
{
  chi_log.LogEvent(source_event_tag, ChiLog::EventType::EVENT_BEGIN);

  const bool apply_wgs_fission_src = (source_flags & APPLY_WGS_FISSION_SOURCE);
  const bool apply_ags_fission_src = (source_flags & APPLY_AGS_FISSION_SOURCE);

  //================================================== Get group setup
  int gs_i = groupset.groups[0].id;
  int gs_f = groupset.groups.back().id;

  int first_grp = groups.front().id;
  int last_grp = groups.back().id;

  const auto& m_to_ell_em_map = groupset.quadrature->GetMomentToHarmonicsIndexMap();

  std::vector<double> default_zero_src(groups.size(), 0.0);

  //================================================== Loop over local cells
  for (auto& cell : grid->local_cells)
  {
    auto& full_cell_view = cell_transport_views[cell.local_id];

    //==================== Obtain xs
    int cell_matid = cell.material_id;
    int xs_id = matid_to_xs_map[cell_matid];

    if ((xs_id < 0) || (xs_id >= material_xs.size()))
    {
      chi_log.Log(LOG_ALLERROR)
          << "Cross-section lookup error\n";
      exit(EXIT_FAILURE);
    }

    auto xs = material_xs[xs_id];

    //======================================== Loop over nodes
    int num_nodes = full_cell_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      size_t ir = full_cell_view.MapDOF(i, 0, 0);
      double* q_mom = &destination_q[ir];
      double* phi_oldp = &phi_old_local[ir];

      //======================================== Loop over groupset groups
      for (int g = gs_i; g <= gs_f; ++g)
      {
        double infission_g = 0.0;
        const bool fission_avail = (xs->is_fissile);

        //=================================== Apply across-groupset fission
        if (fission_avail and apply_ags_fission_src)
        {
          //============================== Loop over groups
          for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
          {
            if (gprime < gs_i or gprime > gs_f)
            {
              // without delayed neutrons
              if (not options.use_precursors)
                infission_g += xs->chi[g] *
                               xs->nu_sigma_f[gprime] *
                               phi_oldp[gprime] / k_eff;

              // with delayed neutrons
              else
              {
                //============================== Prompt fission
                infission_g += xs->chi[g] *
                               xs->nu_prompt_sigma_f[gprime] *
                               phi_oldp[gprime] / k_eff;

                //============================== Delayed fission
                for (size_t j = 0; j < xs->num_precursors; ++j)
                  infission_g += xs->chi_delayed[g][j] *
                                 xs->precursor_yield[j] *
                                 xs->nu_delayed_sigma_f[gprime] *
                                 phi_oldp[gprime] / k_eff;
              }//if use precursors
            }//if across groupset
          }//for gprime
        }//if fission avial


        //=================================== Apply within-groupset fission
        if (fission_avail and apply_wgs_fission_src)
        {
          //============================== Loop over groups
          for (size_t gprime = first_grp; gprime <= last_grp; ++gprime)
          {
            if (gprime >= gs_i and gprime <= gs_f)
            {
              // without delayed neutrons
              if (not options.use_precursors)
                infission_g += xs->chi[g] *
                               xs->nu_sigma_f[gprime] *
                               phi_oldp[gprime] / k_eff;

                // with delayed neutrons
              else
              {
                //========== Prompt fission
                infission_g += xs->chi[g] *
                               xs->nu_prompt_sigma_f[gprime] *
                               phi_oldp[gprime] / k_eff;

                //========== Delayed fission
                for (size_t j = 0; j < xs->num_precursors; ++j)
                  infission_g += xs->chi_delayed[g][j] *
                                 xs->precursor_yield[j] *
                                 xs->nu_delayed_sigma_f[gprime] *
                                 phi_oldp[gprime] / k_eff;
              }//if use precursors
            }//if across groupset
          }//for gprime
        }//if fission avial
        q_mom[g] += infission_g;

      }//for g
    }//for i
  }//for cell

  chi_log.LogEvent(source_event_tag, ChiLog::EventType::EVENT_END);
}