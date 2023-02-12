#include "../lbkes_k_eigenvalue_solver.h"

#include "B_LBSSteadyState/IterativeOperations/sweep_wgs_context.h"
#include "B_LBSSteadyState/IterativeMethods/wgs_linear_solver.h"

#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"
namespace sweep_namespace = chi_mesh::sweep_management;
typedef sweep_namespace::SweepScheduler MainSweepScheduler;
typedef sweep_namespace::SchedulingAlgorithm SchedulingAlgorithm;

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiTimer/chi_timer.h"

#include <iomanip>

using namespace lbs;

//###################################################################
/**Power iterative scheme for k-eigenvalue calculations.
\author Zachary Hardy
*/
void KEigenvalueSolver::PowerIteration()
{
  chi::log.Log()
      << "\n********** Solving k-eigenvalue problem with "
      << "the Power Method.\n";

  phi_old_local_.assign(phi_old_local_.size(), 1.0);

  double F_prev = 1.0;
  k_eff = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  //================================================== Initialize groupsets_
//  for (auto& groupset : groupsets_)
//  {
//    ComputeSweepOrderings(groupset);
//    InitFluxDataStructures(groupset);
//
//    InitWGDSA(groupset);
//    InitTGDSA(groupset);
//  }

  //================================================== Start power iterations
  int nit = 0;
  bool converged = false;
  while (nit < max_iterations)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    // Divide phi_old by k_eff (phi_old gives better init-quess for GMRES)
    for (auto& phi : phi_old_local_) phi /= k_eff;

    //============================================= Loop over groupsets_
    for (auto& groupset : groupsets_)
    {
      //======================================== Precompute the fission source
      q_moments_local_.assign(q_moments_local_.size(), 0.0);
      SetSource(groupset, q_moments_local_,
                PhiOldLocal(),
                APPLY_AGS_FISSION_SOURCES |
                APPLY_WGS_FISSION_SOURCES);

      //======================================== Converge the scattering source
      //                                         with a fixed fission source
      auto sweep_chunk = SetSweepChunk(groupset);

      auto sweep_wgs_context_ptr =
        std::make_shared<SweepWGSContext<Mat, Vec, KSP>>(
          *this, groupset,
          active_set_source_function_,
          APPLY_WGS_SCATTER_SOURCES ,  //lhs_scope
          APPLY_AGS_SCATTER_SOURCES,   //rhs_scope
          true/*with_delayed_psi*/,
          options_.verbose_inner_iterations,
          sweep_chunk);

      WGSLinearSolver<Mat,Vec,KSP> solver(sweep_wgs_context_ptr);
      solver.Setup();
      solver.Solve();

      MPI_Barrier(MPI_COMM_WORLD);
    }//for groupset

    //======================================== Recompute k-eigenvalue
    double F_new = ComputeFissionProduction(phi_new_local_);
    k_eff = F_new / F_prev * k_eff;
    double reactivity = (k_eff - 1.0) / k_eff;

    //======================================== Check convergence, bookkeeping
    k_eff_change = fabs(k_eff - k_eff_prev) / k_eff;
    k_eff_prev = k_eff;
    F_prev = F_new;
    nit += 1;

    if (k_eff_change < std::max(tolerance, 1.0e-12))
      converged = true;

    //======================================== Print iteration summary
    if (options_.verbose_outer_iterations)
    {
      std::stringstream k_iter_info;
      k_iter_info
          << chi::program_timer.GetTimeString() << " "
          << "  Iteration " << std::setw(5) << nit
          << "  k_eff " << std::setw(11) << std::setprecision(7) << k_eff
          << "  k_eff change " << std::setw(12) << k_eff_change
          << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged) k_iter_info << " CONVERGED\n";

      chi::log.Log() << k_iter_info.str();
    }

    if (converged) break;
  }//for k iterations

  //================================================== Cleanup groupsets_
//  for (auto& groupset : groupsets_)
//  {
//    CleanUpWGDSA(groupset);
//    CleanUpTGDSA(groupset);
//
//    ResetSweepOrderings(groupset);
//  }

  //================================================== Print summary
  chi::log.Log() << "\n";
  chi::log.Log()
      << "        Final k-eigenvalue    :        "
      << std::setprecision(7) << k_eff;
  chi::log.Log()
      << "        Final change          :        "
      << std::setprecision(6) << k_eff_change;
  chi::log.Log() << "\n";
}
