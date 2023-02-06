#ifndef CHITECH_GS_CONTEXT_H
#define CHITECH_GS_CONTEXT_H

#include "ChiMath/LinearSolver/linear_solver_context.h"

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

template<class MatType, class VecType>
struct GSContext : public chi_math::LinearSolverContext<MatType, VecType>
{
  LBSGroupset& groupset_;
  SteadyStateSolver& lbs_solver_;
  const SetSourceFunction& set_source_function_;
  const int lhs_src_scope_;
  const int rhs_src_scope_;

  GSContext(LBSGroupset& groupset,
            SteadyStateSolver& lbs_solver,
            const SetSourceFunction& set_source_function,
            int lhs_scope, int rhs_scope) :
            groupset_(groupset),
            lbs_solver_(lbs_solver),
            set_source_function_(set_source_function),
            lhs_src_scope_(lhs_scope),
            rhs_src_scope_(rhs_scope)
  {}

  int MatrixAction(MatType& matrix,
                   VecType& action_vector,
                   VecType& action) override;

  /**This operation applies the inverse of the transform operator in the form
   * Ay = x where the the vector x's underlying implementing is always LBS's
   * q_moments_local vextor.*/
  virtual void ApplyInverseTransportOperator();
};

}//namespace lbs

#endif //CHITECH_GS_CONTEXT_H
