#ifndef LBS_KSP_DATA_CONTEXT_H
#define LBS_KSP_DATA_CONTEXT_H

#include "../lbs_linear_boltzmann_solver.h"

namespace LinearBoltzmann
{

//###################################################################
/**This is a simple data structure of basically pointers to
 * objects needed in matrix free operations.*/
struct KSPDataContext
{
  LinearBoltzmann::Solver& solver;
  LBSGroupset&             groupset;
  Vec&                     operating_vector;
  chi_mesh::sweep_management::SweepScheduler& sweep_scheduler;
  SourceFlags    lhs_scope;
  int64_t last_iteration = -1;

  KSPDataContext(LinearBoltzmann::Solver& in_solver,
                 LBSGroupset& in_groupset,
                 Vec& in_operating_vector,
                 chi_mesh::sweep_management::SweepScheduler& in_sweep_scheduler,
                 SourceFlags in_lhs_scope) :
    solver(in_solver),
    groupset(in_groupset),
    operating_vector(in_operating_vector),
    sweep_scheduler(in_sweep_scheduler),
    lhs_scope(in_lhs_scope) {}
};

}

#endif //LBS_KSP_DATA_CONTEXT_H