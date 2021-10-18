#include "../lbkes_k_eigenvalue_solver.h"

#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

#include <iomanip>

namespace sweep_namespace = chi_mesh::sweep_management;
typedef sweep_namespace::SweepScheduler MainSweepScheduler;
typedef sweep_namespace::SchedulingAlgorithm SchedulingAlgorithm;

using namespace LinearBoltzmann;


//###################################################################
/**Power iterative scheme for k-eigenvalue calculations.
 * Note that this routine currently only works when the problem
 * is defined by a single groupset.
*/
void KEigenvalueSolver::PowerIteration()
{
  chi_log.Log(LOG_0)
      << "\n********** Solving k-eigenvalue problem with "
      << "the Power Method.\n";

  phi_old_local.assign(phi_old_local.size(), 1.0);

  double F_prev = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  //================================================== Start power iterations
  int nit = 0;
  bool converged = false;
  while (nit < max_iterations)
  {
    //============================================= Loop over groupsets
    MPI_Barrier(MPI_COMM_WORLD);
    for (auto& groupset : groupsets)
    {
      ComputeSweepOrderings(groupset);
      InitFluxDataStructures(groupset);

      InitWGDSA(groupset);
      InitTGDSA(groupset);

      //======================================== Setup sweep chunk
      auto sweep_chunk = SetSweepChunk(groupset);
      MainSweepScheduler sweep_scheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
                                         groupset.angle_agg,
                                         *sweep_chunk);

      //======================================== Precompute the fission source
      q_moments_local.assign(q_moments_local.size(), 0.0);
      SetSource(groupset, q_moments_local,
                APPLY_AGS_FISSION_SOURCE |
                APPLY_WGS_FISSION_SOURCE);

      //normalize q by k_eff
      for (auto& q : q_moments_local) q /= k_eff;

      //======================================== Converge the scattering source
      //                                         with a fixed fission source
      if (groupset.iterative_method == IterativeMethod::CLASSICRICHARDSON)
      {
        ClassicRichardson(groupset, sweep_scheduler,
                          APPLY_WGS_SCATTER_SOURCE |
                          APPLY_AGS_SCATTER_SOURCE,
                          options.verbose_inner_iterations);
      }
      else if (groupset.iterative_method == IterativeMethod::GMRES)
      {
        GMRES(groupset, sweep_scheduler,
              APPLY_WGS_SCATTER_SOURCE,
              APPLY_AGS_SCATTER_SOURCE,
              options.verbose_inner_iterations);
      }

      CleanUpWGDSA(groupset);
      CleanUpTGDSA(groupset);

      ResetSweepOrderings(groupset);

      MPI_Barrier(MPI_COMM_WORLD);
    }//for groupset

    //======================================== Recompute k-eigenvalue
    double F_new = ComputeFissionProduction();
    k_eff = F_new / F_prev * k_eff;
    double reactivity = (k_eff - 1.0) / k_eff;

    //======================================== Check convergence, book-keeping
    k_eff_change = fabs(k_eff - k_eff_prev) / k_eff;
    k_eff_prev = k_eff;
    F_prev = F_new;
    nit += 1;

    if (k_eff_change < std::max(tolerance, 1.0e-12))
      converged = true;

    //======================================== Print iteration summary
    if (options.verbose_outer_iterations)
    {
      std::stringstream k_iter_info;
      k_iter_info
          << chi_program_timer.GetTimeString() << " "
          << "  Iteration " << std::setw(5) << nit
          << "  k_eff " << std::setw(10) << k_eff
          << "  k_eff change " << std::setw(10) << k_eff_change
          << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged) k_iter_info << " CONVERGED\n";

      chi_log.Log(LOG_0) << k_iter_info.str();
    }

    if (converged) break;
  }//for k iterations

  //================================================== Print summary
  chi_log.Log(LOG_0) << "\n";
  chi_log.Log(LOG_0)
      << "        Final k-eigenvalue    :        "
      << std::setprecision(6) << k_eff;
  chi_log.Log(LOG_0)
      << "        Final change          :        "
      << std::setprecision(6) << k_eff_change;
  chi_log.Log(LOG_0) << "\n";
}
