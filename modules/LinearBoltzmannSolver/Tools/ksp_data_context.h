#ifndef LBS_KSP_DATA_CONTEXT_H
#define LBS_KSP_DATA_CONTEXT_H

#include "../lbs_linear_boltzmann_solver.h"
#include "LinearBoltzmannSolver/Groupset/lbs_groupset.h"

namespace lbs
{

//###################################################################
/**This is a simple data structure of basically pointers to
 * objects needed in matrix free operations.*/
struct KSPDataContext
{
  lbs::SteadySolver& solver;
  LBSGroupset&             groupset;
  Vec&                     operating_vector;
  chi_mesh::sweep_management::SweepScheduler& sweep_scheduler;
  SourceFlags    lhs_scope;
  int64_t last_iteration = -1;
  double rhs_preconditioned_norm = 0.0;

  const lbs::SetSourceFunction& set_source_function;

  KSPDataContext(lbs::SteadySolver& in_solver,
                 LBSGroupset& in_groupset,
                 Vec& in_operating_vector,
                 chi_mesh::sweep_management::SweepScheduler& in_sweep_scheduler,
                 SourceFlags in_lhs_scope,
                 const lbs::SetSourceFunction& in_set_source_function) :
    solver(in_solver),
    groupset(in_groupset),
    operating_vector(in_operating_vector),
    sweep_scheduler(in_sweep_scheduler),
    lhs_scope(in_lhs_scope),
    set_source_function(in_set_source_function){}
};

}

#endif //LBS_KSP_DATA_CONTEXT_H