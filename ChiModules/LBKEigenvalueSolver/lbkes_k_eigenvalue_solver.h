#ifndef LBKES_K_EIGENVALUE_SOLVER_H
#define LBKES_K_EIGENVALUE_SOLVER_H

#include "LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <string>

namespace lbs
{

/**A k-eigenvalue linear boltzmann transport solver.*/
class KEigenvalueSolver : public lbs::SteadySolver
{
public:
  /**The current k-eigenvalue estimate.*/
  double k_eff = 1.0;

  /**Iterative parameters.*/
  size_t max_iterations = 1000;
  double tolerance = 1.0e-8;

public:
  explicit KEigenvalueSolver(const std::string& in_text_name) :
    lbs::SteadySolver(in_text_name) {}

  void Execute() override;

  // IterativeMethods
  void PowerIteration();

  // Iterative operations
  double ComputeFissionProduction();
};

}

#endif //LBKES_K_EIGENVALUE_SOLVER_H