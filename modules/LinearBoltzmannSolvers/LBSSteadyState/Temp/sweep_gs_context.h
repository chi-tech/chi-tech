#ifndef CHITECH_SWEEP_GS_CONTEXT_H
#define CHITECH_SWEEP_GS_CONTEXT_H


#include "LinearBoltzmannSolvers/LBSSteadyState/Temp/gs_context.h"

#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"
#include "LinearBoltzmannSolvers/LBSSteadyState/Groupset/lbs_groupset.h"

namespace lbs
{

template<class MatType, class VecType, class SolverType>
struct SweepGSContext : public GSContext<MatType,VecType,SolverType>
{
  chi_mesh::sweep_management::SweepScheduler sweep_scheduler_;

  SweepGSContext(LBSGroupset& groupset,
                 SteadyStateSolver& lbs_solver,
                 const SetSourceFunction& set_source_function,
                 int lhs_scope, int rhs_scope,
                 bool with_delayed_psi,
                 bool log_info,
                 std::shared_ptr<chi_mesh::sweep_management::SweepChunk>&
                   sweep_chunk) :
    GSContext<MatType, VecType, SolverType>(groupset,
                                            lbs_solver, set_source_function,
                                            lhs_scope, rhs_scope,
                                            with_delayed_psi,
                                            log_info),
    sweep_scheduler_(
      chi_mesh::sweep_management::SchedulingAlgorithm::DEPTH_OF_GRAPH,
      groupset.angle_agg,
      *sweep_chunk)
  {}

  void PreSetupCallback() override;

  void SetPreconditioner(SolverType& solver) override;

  std::pair<int64_t, int64_t> SystemSize() override;

  void ApplyInverseTransportOperator(int scope) override;

  void PostSolveCallback() override;
};

}//namespace lbs

#endif //CHITECH_SWEEP_GS_CONTEXT_H
