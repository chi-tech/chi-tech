#ifndef LBS_KSP_DATA_CONTEXT_H
#define LBS_KSP_DATA_CONTEXT_H

#include "LBSSteadyState/lbs_linear_boltzmann_solver.h"
#include "LBSSteadyState/Groupset/lbs_groupset.h"

namespace lbs
{

//###################################################################
/**This is a simple data structure of basically pointers to
 * objects needed in matrix free operations.*/
struct KSPDataContext
{
  lbs::SteadyStateSolver& solver;
  LBSGroupset&             groupset;
  chi_mesh::sweep_management::SweepScheduler& sweep_scheduler;
  SourceFlags    lhs_scope;
  int64_t last_iteration = -1;
  double rhs_preconditioned_norm = 0.0;

  const lbs::SetSourceFunction& set_source_function;
  std::vector<double>& phi_old_local;
  std::vector<double>& q_moments_local;
  std::vector<double>& phi_new_local;

  KSPDataContext(lbs::SteadyStateSolver& in_solver,
                 LBSGroupset& in_groupset,
                 chi_mesh::sweep_management::SweepScheduler& in_sweep_scheduler,
                 SourceFlags in_lhs_scope,
                 const lbs::SetSourceFunction& in_set_source_function,
                 std::vector<double>& in_phi_old_local,
                 std::vector<double>& in_q_moments_local,
                 std::vector<double>& in_phi_new_local) :
    solver(in_solver),
    groupset(in_groupset),
    sweep_scheduler(in_sweep_scheduler),
    lhs_scope(in_lhs_scope),
    set_source_function(in_set_source_function),
    phi_old_local(in_phi_old_local),
    q_moments_local(in_q_moments_local),
    phi_new_local(in_phi_new_local)
    {}
};

}

#endif //LBS_KSP_DATA_CONTEXT_H