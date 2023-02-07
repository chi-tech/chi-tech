#ifndef CHITECH_MIP_WGS_CONTEXT_H
#define CHITECH_MIP_WGS_CONTEXT_H

#include "LinearBoltzmannSolvers/LBSSteadyState/Tools/wgs_context.h"

namespace lbs
{

template<class MatType, class VecType, class SolverType>
struct MIPWGSContext : public WGSContext<MatType,VecType,SolverType>
{
  MIPWGSContext(SteadyStateSolver& lbs_solver,
                LBSGroupset& groupset,
                const SetSourceFunction& set_source_function,
                int lhs_scope, int rhs_scope,
                bool log_info) :
    WGSContext<MatType, VecType, SolverType>(lbs_solver,
                                             groupset,
                                             set_source_function,
                                             lhs_scope, rhs_scope,
                                             false,
                                             log_info)
  {}

  void PreSetupCallback() override;

  void SetPreconditioner(SolverType& solver) override;

  std::pair<int64_t, int64_t> SystemSize() override;

  void ApplyInverseTransportOperator(int scope) override;

  void PostSolveCallback() override;
};
}//namespace lbs

#endif //CHITECH_MIP_WGS_CONTEXT_H
