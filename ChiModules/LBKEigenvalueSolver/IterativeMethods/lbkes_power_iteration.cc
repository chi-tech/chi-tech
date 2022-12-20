#include "../lbkes_k_eigenvalue_solver.h"

#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiTimer/chi_timer.h"

#include <iomanip>

namespace sweep_namespace = chi_mesh::sweep_management;
typedef sweep_namespace::SweepScheduler MainSweepScheduler;
typedef sweep_namespace::SchedulingAlgorithm SchedulingAlgorithm;

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

  phi_old_local.assign(phi_old_local.size(), 1.0);

  double F_prev = 1.0;
  k_eff = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  //================================================== Initialize groupsets
  for (auto& groupset : groupsets)
  {
    ComputeSweepOrderings(groupset);
    InitFluxDataStructures(groupset);

    InitWGDSA(groupset);
    InitTGDSA(groupset);
  }

  //================================================== Start power iterations
  int nit = 0;
  bool converged = false;
  while (nit < max_iterations)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    // Divide phi_old by k_eff (phi_old gives better init-quess for GMRES)
    for (auto& phi : phi_old_local) phi /= k_eff;

    //============================================= Loop over groupsets
    for (auto& groupset : groupsets)
    {
      //======================================== Setup sweep chunk
      auto sweep_chunk_ptr = SetSweepChunk(groupset);
      MainSweepScheduler sweep_scheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
                                         groupset.angle_agg,
                                         *sweep_chunk_ptr);

      //======================================== Precompute the fission source
      q_moments_local.assign(q_moments_local.size(), 0.0);
      SetSource(groupset, q_moments_local,
                APPLY_AGS_FISSION_SOURCE |
                APPLY_WGS_FISSION_SOURCE);

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
      else if (groupset.iterative_method == IterativeMethod::KRYLOV_RICHARDSON or
               groupset.iterative_method == IterativeMethod::KRYLOV_GMRES or
               groupset.iterative_method == IterativeMethod::KRYLOV_BICGSTAB)
      {
        Krylov(groupset, sweep_scheduler,
               APPLY_WGS_SCATTER_SOURCE,
               APPLY_AGS_SCATTER_SOURCE,
               options.verbose_inner_iterations);
      }

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
          << chi::program_timer.GetTimeString() << " "
          << "  Iteration " << std::setw(5) << nit
          << "  k_eff " << std::setw(10) << k_eff
          << "  k_eff change " << std::setw(10) << k_eff_change
          << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged) k_iter_info << " CONVERGED\n";

      chi::log.Log() << k_iter_info.str();
    }

    if (converged) break;
  }//for k iterations

  //================================================== Cleanup groupsets
  for (auto& groupset : groupsets)
  {
    CleanUpWGDSA(groupset);
    CleanUpTGDSA(groupset);

    ResetSweepOrderings(groupset);
  }

  //================================================== Print summary
  chi::log.Log() << "\n";
  chi::log.Log()
      << "        Final k-eigenvalue    :        "
      << std::setprecision(6) << k_eff;
  chi::log.Log()
      << "        Final change          :        "
      << std::setprecision(6) << k_eff_change;
  chi::log.Log() << "\n";
}
