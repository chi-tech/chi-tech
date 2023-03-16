#include "lbkes_k_eigenvalue_solver.h"
#include "A_LBSSolver/IterativeMethods/poweriteration_keigen.h"
#include "A_LBSSolver/IterativeMethods/nl_keigen_ags_solver.h"

#include "chi_log.h"

#include <petscsnes.h>

//###################################################################
/**Execute a k-eigenvalue linear boltzmann solver.*/
void lbs::DiscOrdKEigenvalueSolver::Execute()
{
  const std::string fname = "lbs::DiscOrdKEigenvalueSolver::Execute";

  //======================================== Solve the k-eigenvalue problem
  if (k_eigen_method_ == "power")
  {
    chi::log.Log()
      << "\n********** Solving k-eigenvalue problem with "
      << "the Power Iteration method.\n";

    lbs::PowerIterationKEigen(*this, tolerance_, max_iterations_, k_eff_);
  }
  else if (k_eigen_method_ == "nonlinear")
  {
    chi::log.Log()
      << "\n********** Solving k-eigenvalue problem with "
      << "the non-linear k-eigenvalue method.\n";

    SetPhiVectorScalarValues(PhiOldLocal(), 1.0);

    auto context = std::make_shared<NLKEigenAGSContext<Vec,SNES>>(*this);

    NLKEigenvalueAGSSolver<Mat,Vec,SNES> nl_solver(SNESNEWTONLS, context);

    auto& tolerances = nl_solver.ToleranceOptions();

    tolerances.nl_absolute_tol = tolerance_;
    tolerances.nl_solution_tol = tolerance_;
    tolerances.nl_max_iterations = max_iterations_;

    const auto& front_gs = groupsets_.front();
    tolerances.l_absolute_tol = front_gs.residual_tolerance_;
    tolerances.l_relative_tol = front_gs.residual_tolerance_;
    tolerances.l_max_iterations = front_gs.max_iterations_;
    tolerances.l_gmres_restart_interval = front_gs.gmres_restart_intvl_;

    nl_solver.Setup();
    nl_solver.Solve();
  }
  else
    throw std::invalid_argument(fname + ": Unsupported k_eigen_method. "
      "Must be any of \"power\", \"nonlinear\"");

  //======================================== Initialize the precursors
  if (options_.use_precursors)
  {
    ComputePrecursors();
    for (auto& v : precursor_new_local_) v /= k_eff_;
  }

  UpdateFieldFunctions();

  chi::log.Log()
      << "LinearBoltzmann::KEigenvalueSolver execution completed\n\n";
}