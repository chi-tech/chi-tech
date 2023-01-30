#include "lbsadj_solver.h"

//###################################################################
/**Constructor.*/
lbs::SteadyStateAdjointSolver::SteadyStateAdjointSolver(const std::string &solver_name) :
  lbs::SteadyStateSolver(solver_name)
{
  basic_options.AddOption<std::string>("REFERENCE_RF", std::string());
}