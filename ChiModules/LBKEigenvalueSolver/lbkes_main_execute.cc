#include "lbkes_k_eigenvalue_solver.h"

#include <chi_log.h>
;

#include <iomanip>

using namespace lbs;

//###################################################################
/**Execute a k-eigenvalue linear boltzmann solver.*/
void KEigenvalueSolver::Execute()
{
  //======================================== Solve the k-eigenvalue problem
  PowerIteration();

  //======================================== Initialize the precursors
  if (options.use_precursors)
  {
    ComputePrecursors();
    for (auto& v : precursor_new_local) v /= k_eff;
  }

  chi::log.Log()
      << "LinearBoltzmann::KEigenvalueSolver execution completed\n\n";
}