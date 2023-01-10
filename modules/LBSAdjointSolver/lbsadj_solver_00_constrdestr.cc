#include "lbsadj_solver.h"

//###################################################################
/**Constructor.*/
lbs_adjoint::AdjointSolver::AdjointSolver(const std::string &solver_name) :
  lbs::SteadySolver(solver_name)
{
  basic_options.AddOption<std::string>("REFERENCE_RF", std::string());
}