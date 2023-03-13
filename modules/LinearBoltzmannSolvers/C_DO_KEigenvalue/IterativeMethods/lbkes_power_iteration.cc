#include "../lbkes_k_eigenvalue_solver.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiTimer/chi_timer.h"

#include <iomanip>

using namespace lbs;

//###################################################################
/**Power iterative scheme for k-eigenvalue calculations.
\author Zachary Hardy
*/
void DiscOrdKEigenvalueSolver::PowerIteration()
{
  chi::log.Log()
      << "\n********** Solving k-eigenvalue problem with "
      << "the Power Method.\n";

  SetPhiVectorScalarValues(phi_old_local_, 1.0);

  double F_prev = 1.0;
  k_eff_ = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  //================================================== Start power iterations
  primary_ags_solver_->SetVerbosity(options_.verbose_ags_iterations);
  int nit = 0;
  bool converged = false;
  while (nit < max_iterations_)
  {
    // Divide phi_old by k_eff (phi_old gives better init-quess for GMRES)
    for (auto& phi : phi_old_local_) phi /= k_eff_;

    q_moments_local_.assign(q_moments_local_.size(), 0.0);
    for (auto& groupset : groupsets_)
    {
      active_set_source_function_(groupset, q_moments_local_,
                                  PhiOldLocal(),
                                  APPLY_AGS_FISSION_SOURCES |
                                  APPLY_WGS_FISSION_SOURCES);
      groupset.angle_agg_.ZeroIncomingDelayedPsi();
    }


    primary_ags_solver_->Setup();
    primary_ags_solver_->Solve();

    //======================================== Recompute k-eigenvalue
    double F_new = ComputeFissionProduction(phi_new_local_);
    k_eff_ = F_new / F_prev * k_eff_;
    double reactivity = (k_eff_ - 1.0) / k_eff_;

    //======================================== Check convergence, bookkeeping
    k_eff_change = fabs(k_eff_ - k_eff_prev) / k_eff_;
    k_eff_prev = k_eff_;
    F_prev = F_new;
    nit += 1;

    if (k_eff_change < std::max(tolerance_, 1.0e-12))
      converged = true;

    //======================================== Print iteration summary
    if (options_.verbose_outer_iterations)
    {
      std::stringstream k_iter_info;
      k_iter_info
        << chi::program_timer.GetTimeString() << " "
        << "  Iteration " << std::setw(5) << nit
        << "  k_eff " << std::setw(11) << std::setprecision(7) << k_eff_
          << "  k_eff change " << std::setw(12) << k_eff_change
          << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged) k_iter_info << " CONVERGED\n";

      chi::log.Log() << k_iter_info.str();
    }

    if (converged) break;
  }//for k iterations

  //================================================== Print summary
  chi::log.Log() << "\n";
  chi::log.Log()
    << "        Final k-eigenvalue    :        "
    << std::setprecision(7) << k_eff_;
  chi::log.Log()
      << "        Final change          :        "
      << std::setprecision(6) << k_eff_change;
  chi::log.Log() << "\n";
}
