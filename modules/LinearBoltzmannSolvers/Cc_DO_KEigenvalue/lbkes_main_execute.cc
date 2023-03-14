#include "lbkes_k_eigenvalue_solver.h"
#include "A_LBSSolver/IterativeMethods/poweriteration_keigen.h"
#include "A_LBSSolver/IterativeMethods/nonlinear_keigen.h"

#include "chi_log.h"

//###################################################################
/**Execute a k-eigenvalue linear boltzmann solver.*/
void lbs::DiscOrdKEigenvalueSolver::Execute()
{
  const std::string fname = "lbs::DiscOrdKEigenvalueSolver::Execute";

  //======================================== Solve the k-eigenvalue problem
  if (k_eigen_method_ == "power")
  {
    lbs::PowerIterationKEigen(*this, tolerance_, max_iterations_, k_eff_);
  }
  else if (k_eigen_method_ == "nonlinear")
  {
    lbs::NonLinearKEigen(*this, tolerance_, max_iterations_, k_eff_);
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