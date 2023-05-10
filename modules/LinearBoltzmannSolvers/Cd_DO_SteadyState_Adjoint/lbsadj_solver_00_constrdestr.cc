#include "lbsadj_solver.h"

//###################################################################
/**Constructor.*/
lbs::DiscOrdSteadyStateAdjointSolver::
  DiscOrdSteadyStateAdjointSolver(const std::string &solver_name) :
    lbs::DiscreteOrdinatesSolver(solver_name)
{
  basic_options_.AddOption<std::string>("REFERENCE_RF", std::string());
}

/**Returns the list of volumetric response functions.*/
const std::vector<lbs::DiscOrdSteadyStateAdjointSolver::RespFuncAndSubs>&
  lbs::DiscOrdSteadyStateAdjointSolver::GetResponseFunctions() const
{
  return response_functions_;
}
