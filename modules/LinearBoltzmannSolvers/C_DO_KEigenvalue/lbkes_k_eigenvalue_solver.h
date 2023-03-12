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
  double k_eff_ = 1.0;

  /**Iterative parameters.*/
  size_t max_iterations_ = 1000;
  double tolerance_ = 1.0e-8;
  std::string k_eigen_method_ = "power";

public:
  DiscOrdKEigenvalueSolver (const DiscOrdKEigenvalueSolver&) = delete;
  DiscOrdKEigenvalueSolver& operator= (const DiscOrdKEigenvalueSolver&) = delete;

  explicit DiscOrdKEigenvalueSolver(const std::string& in_text_name) :
    lbs::DiscOrdSteadyStateSolver(in_text_name) {}

  double GetKeff() const {return k_eff_;}
  void SetKeff(double k_eff) {k_eff_ = k_eff;}

protected:
  void InitializeWGSSolvers() override;

public:
  void Execute() override;

  // IterativeMethods
  void PowerIteration();
  int NonLinearKEigen();


};

}

#endif //LBKES_K_EIGENVALUE_SOLVER_H