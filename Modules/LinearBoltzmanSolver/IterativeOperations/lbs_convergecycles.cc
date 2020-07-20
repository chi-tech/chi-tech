#include "../lbs_linear_boltzman_solver.h"
#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include <ChiTimer/chi_timer.h>
#include <ChiLog/chi_log.h>

extern ChiLog& chi_log;
extern ChiTimer chi_program_timer;

#include <iomanip>

//###################################################################
/**Performs iterations without source updates to converge cyclic
 * dependencies.*/
void LinearBoltzman::Solver::ConvergeCycles(
  MainSweepScheduler& sweepScheduler,
  SweepChunk* sweep_chunk,
  LBSGroupset *groupset,
  bool convergence_opp_refl_bndries,
  bool apply_latest_convergence_metric,
  double cyclic_tolerance,
  size_t cyclic_max_iter)
{
  auto REFLECTING_BNDRY = chi_mesh::sweep_management::BoundaryType::REFLECTING;
  typedef chi_mesh::sweep_management::BoundaryReflecting TBndryReflecting;

  double convergence_metric = cyclic_tolerance;
  if (apply_latest_convergence_metric)
    convergence_metric =
      std::max(cyclic_tolerance,0.8*groupset->latest_convergence_metric);

  double max_pw_change = groupset->angle_agg->GetDelayedPsiNorm();
  if (convergence_opp_refl_bndries)
    for (auto bndry : sweep_boundaries)
      if (bndry->Type() == REFLECTING_BNDRY)
        max_pw_change += ((TBndryReflecting*)bndry)->pw_change;

  if (max_pw_change<convergence_metric)
    return;

  std::vector<double> temp_phi_old(phi_old_local.size(),0.0);
  DisAssembleVectorLocalToLocal(groupset,phi_old_local.data(), temp_phi_old.data());
  DisAssembleVectorLocalToLocal(groupset,phi_new_local.data(), phi_old_local.data());

  bool cycles_converged = false;
  for (int k=0; k<cyclic_max_iter; k++)
  {
    phi_new_local.assign(phi_new_local.size(),0.0); //Ensure phi_new=0.0
    sweepScheduler.Sweep(sweep_chunk);

    max_pw_change = ComputePiecewiseChange(groupset)*
                    groupset->angle_agg->GetDelayedPsiNorm();

    if (convergence_opp_refl_bndries)
      for (auto bndry : sweep_boundaries)
        if (bndry->Type() == REFLECTING_BNDRY)
          max_pw_change += ((TBndryReflecting*)bndry)->pw_change;

    DisAssembleVectorLocalToLocal(groupset,phi_new_local.data(), phi_old_local.data());

    if (apply_latest_convergence_metric)
      convergence_metric =
        std::max(cyclic_tolerance,0.8*groupset->latest_convergence_metric);

    if (max_pw_change<convergence_metric)
      cycles_converged = true;

    std::stringstream iter_info;
    iter_info
      << chi_program_timer.GetTimeString()
      << " Cyclic iteration " << std::setw(5) << k
      << " Point-wise change " << std::setw(14) << max_pw_change;

    if (cycles_converged)
      iter_info << " CONVERGED\n";

    chi_log.Log(LOG_0) << iter_info.str();

    if (cycles_converged)
      break;
  }

  DisAssembleVectorLocalToLocal(groupset,temp_phi_old.data(), phi_old_local.data());

  groupset->angle_agg->ResetReflectingBCs();
}