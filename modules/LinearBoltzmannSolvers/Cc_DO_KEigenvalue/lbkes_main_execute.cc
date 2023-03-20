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

  auto k_eigen_method = basic_options_("K_EIGEN_METHOD").StringValue();
  bool reset_solution = basic_options_("K_EIGEN_RESET_SOLUTION").BoolValue();

  //======================================== Solve the k-eigenvalue problem
  if (k_eigen_method == "power")
  {
    chi::log.Log()
      << "\n********** Solving k-eigenvalue problem with "
      << "the Power Iteration method.\n";

    auto pi_max_iters =
      static_cast<int>(basic_options_("PI_MAX_ITS").IntegerValue());
    auto pi_k_tolerance =
      basic_options_("PI_K_TOL").FloatValue();

    if (reset_solution) SetPhiVectorScalarValues(PhiOldLocal(), 1.0);

    lbs::PowerIterationKEigen(*this, pi_k_tolerance, pi_max_iters, k_eff_);
  }
  else if (k_eigen_method == "power1")
  {
    chi::log.Log()
      << "\n********** Solving k-eigenvalue problem with "
      << "the Power Iteration method.\n";

    auto pi_max_iters =
      static_cast<int>(basic_options_("PI_MAX_ITS").IntegerValue());
    auto pi_k_tolerance =
      basic_options_("PI_K_TOL").FloatValue();

    if (reset_solution) SetPhiVectorScalarValues(PhiOldLocal(), 1.0);

    lbs::PowerIterationKEigen1(*this, pi_k_tolerance, pi_max_iters, k_eff_);
  }
  else if (k_eigen_method == "power2")
  {
    chi::log.Log()
      << "\n********** Solving k-eigenvalue problem with "
      << "the Power Iteration method.\n";

    auto pi_max_iters =
      static_cast<int>(basic_options_("PI_MAX_ITS").IntegerValue());
    auto pi_k_tolerance =
      basic_options_("PI_K_TOL").FloatValue();

    if (reset_solution) SetPhiVectorScalarValues(PhiOldLocal(), 1.0);

    lbs::PowerIterationKEigen2(*this, pi_k_tolerance, pi_max_iters, k_eff_);
  }
  else if (k_eigen_method == "nonlinear")
  {
    chi::log.Log()
      << "\n********** Solving k-eigenvalue problem with "
      << "the non-linear k-eigenvalue method.\n";

    if (reset_solution) SetPhiVectorScalarValues(PhiOldLocal(), 1.0);

    auto context = std::make_shared<NLKEigenAGSContext<Vec,SNES>>(*this);

    NLKEigenvalueAGSSolver<Mat,Vec,SNES> nl_solver(SNESNEWTONLS, context);

    auto& tolerances = nl_solver.ToleranceOptions();

    tolerances.nl_absolute_tol =
      basic_options_("NLK_ABS_TOL").FloatValue();
    tolerances.nl_solution_tol =
      basic_options_("NLK_SOL_TOL").FloatValue();
    tolerances.nl_max_iterations =
      static_cast<int>(basic_options_("NLK_MAX_ITS").IntegerValue());

    tolerances.l_relative_tol =
      basic_options_("NLK_L_REL_TOL").FloatValue();
    tolerances.l_absolute_tol =
      basic_options_("NLK_L_ABS_TOL").FloatValue();
    tolerances.l_divergence_tol =
      basic_options_("NLK_L_DIV_TOL").FloatValue();
    tolerances.l_max_iterations =
      static_cast<int>(basic_options_("NLK_L_MAX_ITS").IntegerValue());
    tolerances.l_gmres_restart_interval =
      static_cast<int>(basic_options_("NLK_GMRES_RESTART_INTVL").IntegerValue());
    tolerances.l_gmres_breakdown_tol =
      basic_options_("NLK_GMRES_BRKDN_TOL").FloatValue();

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