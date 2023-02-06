#ifndef CHITECH_GSLINSOLVEBASE_H
#define CHITECH_GSLINSOLVEBASE_H

#include "gs_context.h"

#include "ChiMath/LinearSolver/linear_solver.h"

#include <memory>
#include <vector>
#include <functional>

namespace lbs
{
class LBSGroupset;
class SteadyStateSolver;

}

typedef std::function<void(lbs::LBSGroupset&,
                           std::vector<double>&,
                           int)> SetSourceFunction;

namespace lbs
{

template<class MatType, class VecType, class SolverType>
class GSLinearSolver : public chi_math::LinearSolver<MatType, VecType, SolverType>
{
protected:
  SteadyStateSolver& lbs_solver_;
  LBSGroupset& groupset_;
  const SetSourceFunction& set_source_function_;
  const int lhs_src_scope_;
  const int rhs_src_scope_;

  std::vector<double> saved_q_moments_local_;
  bool log_info = true;

public:
  GSLinearSolver(SteadyStateSolver& lbs_solver,
                 LBSGroupset& groupset,
                 const std::string iterative_method,
                 const SetSourceFunction& set_source_function,
                 int lhs_scope, int rhs_scope) :
    chi_math::LinearSolver<MatType, VecType, SolverType>(iterative_method),
    lbs_solver_(lbs_solver),
    groupset_(groupset),
    set_source_function_(set_source_function),
    lhs_src_scope_(lhs_scope),
    rhs_src_scope_(rhs_scope)
  {}
//  virtual void Setup();
  void PreSetupCallback() override;
//  virtual void SetOptions();
  void SetSolverContext() override;
  void SetConvergenceTest() override;
//  virtual void SetMonitor();
//  virtual void SetPreconditioner();

  virtual void SetSystemMatrix();
//  virtual void PostSetupCallback();

//  virtual void Solve();
//  virtual void PreSolveCallback();
  void SetRHS() override;
  void SetInitialGuess() override;
  void PostSolveCallback() override;
};

}//namespace lbs

#endif //CHITECH_GSLINSOLVEBASE_H
