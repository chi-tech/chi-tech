#include "mg_diffusion_solver.h"
#include "chi_runtime.h"
#include "chi_log.h"

//========================================================== Solve 1g problem
void mg_diffusion::Solver::SolveOneGroupProblem(const unsigned int g)
{
  chi::log.Log() << "Solving group: " << g;

  KSPSetOperators(petsc_solver.ksp, A[g], A[g]);

  KSPSolve(petsc_solver.ksp,b,x[g]);
 
  chi::log.Log() << "Done solving group" << g;
}