#include "../lbs_linear_boltzman_solver.h"
#include "ChiMesh/SweepUtilities/chi_sweepscheduler.h"

#include <ChiTimer/chi_timer.h>
#include <ChiLog/chi_log.h>


extern ChiLog chi_log;
extern ChiTimer chi_program_timer;

//###################################################################
/**Performs iterations without source updates to converge cyclic
 * dependencies.*/
void LinearBoltzman::Solver::ConvergeCycles(
  MainSweepScheduler& sweepScheduler,
  SweepChunk* sweep_chunk,
  LBSGroupset *groupset)
{
  double max_pw_change = groupset->angle_agg->GetDelayedPsiNorm();
  if (max_pw_change<std::max(1.0e-8,1.0e-10))
    return;

  bool cycles_converged = false;
  for (int k=0; k<50; k++)
  {
    phi_new_local.assign(phi_new_local.size(),0.0); //Ensure phi_new=0.0
    sweepScheduler.Sweep(sweep_chunk);
    max_pw_change = groupset->angle_agg->GetDelayedPsiNorm();

    if (max_pw_change<std::max(1.0e-8,1.0e-10))
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
}