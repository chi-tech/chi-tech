#ifndef LBKES_K_EIGENVALUE_SOLVER_H
#define LBKES_K_EIGENVALUE_SOLVER_H

#include "B_DO_SteadyState/lbs_DO_steady_state.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <string>

namespace lbs
{

/**A k-eigenvalue solver based on the linear boltzmann transport solver.
\author Zachary Hardy.*/
class DiscOrdKEigenvalueSolver : public lbs::DiscOrdSteadyStateSolver
{
public:
  /**The current k-eigenvalue estimate.*/
  double k_eff = 1.0;

  /**Iterative parameters.*/
  size_t max_iterations = 1000;
  double tolerance = 1.0e-8;

public:
  DiscOrdKEigenvalueSolver (const DiscOrdKEigenvalueSolver&) = delete;
  DiscOrdKEigenvalueSolver& operator= (const DiscOrdKEigenvalueSolver&) = delete;

  explicit DiscOrdKEigenvalueSolver(const std::string& in_text_name) :
    lbs::DiscOrdSteadyStateSolver(in_text_name) {}

  void Execute() override;

  // IterativeMethods
  void PowerIteration();


};

}

#endif //LBKES_K_EIGENVALUE_SOLVER_H