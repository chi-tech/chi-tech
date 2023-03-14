#ifndef LBKES_K_EIGENVALUE_SOLVER_H
#define LBKES_K_EIGENVALUE_SOLVER_H

#include "B_DO_Solver/lbs_discrete_ordinates_solver.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include <string>

namespace lbs
{

/**A k-eigenvalue solver based on the linear boltzmann transport solver.
\author Zachary Hardy.*/
class DiscOrdKEigenvalueSolver : public lbs::LBSDiscreteOrdinatesSolver
{
private:
  /**The current k-eigenvalue estimate.*/
  double k_eff_ = 1.0;

  /**Iterative parameters.*/
  int max_iterations_ = 1000;
  double tolerance_ = 1.0e-8;
  std::string k_eigen_method_ = "power";

public:
  explicit DiscOrdKEigenvalueSolver(const std::string& in_text_name) :
    lbs::LBSDiscreteOrdinatesSolver(in_text_name) {}

  DiscOrdKEigenvalueSolver (const DiscOrdKEigenvalueSolver&) = delete;
  DiscOrdKEigenvalueSolver& operator= (const DiscOrdKEigenvalueSolver&) = delete;

  double GetKeff() const {return k_eff_;}
  void SetKeff(double k_eff) {k_eff_ = k_eff;}

  int GetMaxIterations() const {return max_iterations_;}
  void SetMaxIterations(int max_iterations) {max_iterations_ = max_iterations;}

  double GetTolerance() const {return tolerance_;}
  void SetTolerance(double tolerance) {tolerance_ = tolerance;}

  std::string GetKEigenMethod() const {return k_eigen_method_;}
  void SetKEigenMethod(std::string k_eigen_method)
  {k_eigen_method_ = std::move(k_eigen_method);}

public:
  void Execute() override;
};

}

#endif //LBKES_K_EIGENVALUE_SOLVER_H