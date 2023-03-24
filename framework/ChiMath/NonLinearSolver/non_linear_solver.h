#ifndef CHITECH_CHI_MATH_NON_LINEAR_SOLVER_H
#define CHITECH_CHI_MATH_NON_LINEAR_SOLVER_H

#include "nl_solver_context.h"

#include <string>

#include <memory>
#include <utility>

namespace chi_math
{

/**Implementation of a general non-linear solver.*/
template<class MatType, class VecType, class SolverType>
class NonLinearSolver
{
public:
  typedef NonLinearSolverContext<VecType,SolverType> NLSolverContext;
  typedef std::shared_ptr<NLSolverContext> NLSolverContextPtr;

protected:
  const std::string solver_name_;
  const std::string nl_method_;
  std::string linear_method_ = "gmres";

  NLSolverContextPtr context_ptr_ = nullptr;

  MatType    J_;
  VecType    r_;
  VecType    x_;
  SolverType nl_solver_;

private:
  bool system_set_ = false;

protected:
  int64_t num_local_dofs_ = 0;
  int64_t num_globl_dofs_ = 0;

  struct ToleranceOptions
  {
    double nl_relative_tol = 1.0e-50;
    double nl_absolute_tol = 1.0e-6;
    double nl_solution_tol = 1.0e-6;
    int    nl_max_iterations = 50;
    int    nl_max_r_evaluations = -1;
    int    l_max_failed_iterations = 1000;
    double l_relative_tol = 1.0e-50;
    double l_absolute_tol = 1.0e-6;
    double l_divergence_tol = 1.0e6;
    int    l_max_iterations = 100;
    int    l_gmres_restart_interval = 30;
    double l_gmres_breakdown_tol = 1.0e6;
  }tolerance_options_;

protected:
  bool IsSystemSet() const {return system_set_;}

public:
  NonLinearSolver(const std::string& nl_method,
                  NLSolverContextPtr context_ptr) :
                  solver_name_(nl_method),
                  nl_method_(nl_method),
                  context_ptr_(context_ptr)
  {}
  virtual ~NonLinearSolver();

  ToleranceOptions& ToleranceOptions()
  { return tolerance_options_; }
  void ApplyToleranceOptions();

  NLSolverContextPtr& GetContext()
  { return context_ptr_; }

protected:
  virtual void PreSetupCallback();
  virtual void SetOptions();
  virtual void SetSolverContext();
  virtual void SetConvergenceTest();
  virtual void SetMonitor();
  virtual void SetPreconditioner();

  virtual void SetSystemSize() = 0;
  virtual void SetSystem() = 0;
  virtual void SetFunction() = 0;
  virtual void SetJacobian() = 0;
  virtual void PostSetupCallback();
public:
  virtual void Setup();

protected:
  virtual void PreSolveCallback();
  virtual void SetInitialGuess() = 0;
  virtual void PostSolveCallback();
public:
  virtual void Solve();

};


}//namespace chi_math

#endif //CHITECH_CHI_MATH_NON_LINEAR_SOLVER_H
