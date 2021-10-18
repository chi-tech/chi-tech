#include "lbkes_k_eigenvalue_solver.h"

#include <chi_log.h>
extern ChiLog& chi_log;

#include <iomanip>

using namespace LinearBoltzmann;

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

  chi_log.Log(LOG_0)
      << "LinearBoltzmann::KEigenvalueSolver execution completed\n\n";
}