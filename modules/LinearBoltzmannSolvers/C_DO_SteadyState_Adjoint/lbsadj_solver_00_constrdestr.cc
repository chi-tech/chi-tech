#include "lbsadj_solver.h"

//###################################################################
/**Constructor.*/
lbs::DiscOrdSteadyStateAdjointSolver::DiscOrdSteadyStateAdjointSolver(const std::string &solver_name) :
  lbs::DiscOrdSteadyStateSolver(solver_name)
{
  basic_options_.AddOption<std::string>("REFERENCE_RF", std::string());
}