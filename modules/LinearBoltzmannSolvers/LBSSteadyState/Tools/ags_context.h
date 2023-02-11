#ifndef CHITECH_AGS_CONTEXT_H
#define CHITECH_AGS_CONTEXT_H

#include "ChiMath/LinearSolver/linear_solver_context.h"
#include "ChiMath/LinearSolver/linear_solver.h"

#include <vector>
#include <memory>

namespace lbs
{
  class LBSGroupset;
  class SteadyStateSolver;
}

namespace lbs
{

template<class MatType, class VecType, class SolverType>
struct AGSContext : public chi_math::LinearSolverContext<MatType, VecType>
{
  typedef chi_math::LinearSolver<MatType, VecType, SolverType> LinSolveBaseType;
  typedef std::shared_ptr<LinSolveBaseType> LinSolveBaseTypePtr;
  SteadyStateSolver& lbs_solver_;
  std::vector<LinSolveBaseTypePtr> sub_solvers_list_;

  explicit
  AGSContext(SteadyStateSolver& lbs_solver,
             std::vector<LinSolveBaseTypePtr> sub_solvers_list) :
    lbs_solver_(lbs_solver),
    sub_solvers_list_(std::move(sub_solvers_list))
  {}

  std::pair<int64_t,int64_t> SystemSize();

  virtual void SetPreconditioner(SolverType& solver);

  int MatrixAction(MatType& matrix, VecType& vector, VecType& action) override;

};

}//namespace lbs

#endif //CHITECH_AGS_CONTEXT_H