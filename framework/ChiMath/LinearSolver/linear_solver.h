#ifndef CHITECH_LINEAR_SOLVER_H
#define CHITECH_LINEAR_SOLVER_H

#include "linear_solver_context.h"

#include <string>
#include <utility>

#include <memory>

namespace chi_math
{
template<class MatType, class VecType>
struct LinearSolverContext;

/**Implementation of a general linear solver.*/
template<class MatType, class VecType, class SolverType>
class LinearSolver
{
protected:
  const std::string solver_name_;
  const std::string iterative_method_;

  std::shared_ptr<LinearSolverContext<MatType,VecType>> context_ptr_ = nullptr;

  MatType A_;
  VecType b_;
  VecType x_;
  SolverType solver_;

  int64_t num_local_dofs_ = 0;
  int64_t num_globl_dofs_ = 0;

public:
  struct ToleranceOptions
  {
    double residual_relative_   = 1.0e-50;
    double residual_absolute_   = 1.0e-6;
    double residual_divergence_ = 1.0e6;
    int    maximum_iterations_  = 100;
    int    gmres_restart_interval = 100;
    double gmres_breakdown_tolerance = 1.0e6;
  }tolerance_options_;

public:
  explicit LinearSolver(const std::string& iterative_method) :
    solver_name_(iterative_method),
    iterative_method_(iterative_method) {}

  explicit LinearSolver(std::string  solver_name,
                        std::string  iterative_method) :
    solver_name_(std::move(solver_name)),
    iterative_method_(std::move(iterative_method)) {}

  virtual void Setup();
  virtual void PreSetupCallback();
  virtual void SetOptions();
  virtual void SetSolverContext();
  virtual void SetConvergenceTest();
  virtual void SetMonitor();
  virtual void SetPreconditioner();

  virtual void SetSystemMatrix() = 0; //TODO: Rename to SetSystem
  virtual void PostSetupCallback();

  virtual void Solve();
  virtual void PreSolveCallback();
  virtual void SetRHS() = 0;
  virtual void SetInitialGuess() = 0;
  virtual void PostSolveCallback();
};

}//namespace chi_math

#endif //CHITECH_LINEAR_SOLVER_H
