#ifndef CHITECH_SWEEP_GS_CONTEXT_H
#define CHITECH_SWEEP_GS_CONTEXT_H

#include <petscksp.h>
#include "LinearBoltzmannSolvers/LBSSteadyState/Temp/gs_context.h"

#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

namespace lbs
{

template<class MatType, class VecType>
struct SweepGSContext : public GSContext<MatType,VecType>
{
  chi_mesh::sweep_management::SweepScheduler sweep_scheduler_;

  SweepGSContext(LBSGroupset& groupset,
                 SteadyStateSolver& lbs_solver,
                 const SetSourceFunction& set_source_function,
                 int lhs_scope, int rhs_scope,
                 bool with_delayed_psi,
                 std::shared_ptr<chi_mesh::sweep_management::SweepChunk>&
                   sweep_chunk) :
    GSContext<MatType, VecType>(groupset, lbs_solver, set_source_function,
                                lhs_scope, rhs_scope, with_delayed_psi)
  {}
};

}//namespace lbs

#endif //CHITECH_SWEEP_GS_CONTEXT_H
