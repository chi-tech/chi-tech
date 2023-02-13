#ifndef CHITECH_SWEEP_WGS_CONTEXT_H
#define CHITECH_SWEEP_WGS_CONTEXT_H

#include "A_LBSSolver/Tools/wgs_context.h"

#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

namespace lbs
{
  class DiscOrdSteadyStateSolver;
}

namespace lbs
{

template<class MatType, class VecType, class SolverType>
struct SweepWGSContext : public WGSContext<MatType,VecType,SolverType>
{
  std::shared_ptr<chi_mesh::sweep_management::SweepChunk> sweep_chunk_;
  chi_mesh::sweep_management::SweepScheduler sweep_scheduler_;

  DiscOrdSteadyStateSolver& lbs_ss_solver_;

  SweepWGSContext(DiscOrdSteadyStateSolver& lbs_solver,
                  LBSGroupset& groupset,
                  const SetSourceFunction& set_source_function,
                  int lhs_scope, int rhs_scope,
                  bool with_delayed_psi,
                  bool log_info,
                  std::shared_ptr<chi_mesh::sweep_management::SweepChunk>
                   sweep_chunk) :
    WGSContext<MatType, VecType, SolverType>(lbs_solver,
                                             groupset,
                                             set_source_function,
                                             lhs_scope, rhs_scope,
                                             with_delayed_psi,
                                             log_info),
    sweep_chunk_(std::move(sweep_chunk)),
    sweep_scheduler_(
      chi_mesh::sweep_management::SchedulingAlgorithm::DEPTH_OF_GRAPH,
      groupset.angle_agg,
      *sweep_chunk_),
    lbs_ss_solver_(lbs_solver)
  {}

  void PreSetupCallback() override;

  void SetPreconditioner(SolverType& solver) override;

  std::pair<int64_t, int64_t> SystemSize() override;

  void ApplyInverseTransportOperator(int scope) override;

  void PostSolveCallback() override;
};

}//namespace lbs

#endif //CHITECH_SWEEP_WGS_CONTEXT_H
