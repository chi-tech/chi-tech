#include "diffusion_solver.h"

#include <chi_log.h>
extern ChiLog& chi_log;

chi_diffusion::Solver::Solver()
{}

chi_diffusion::Solver::Solver(std::string in_solver_name):Solver()
{
  solver_name = std::move(in_solver_name);
}

chi_diffusion::Solver::~Solver()
{
  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << "Cleaning up diffusion solver: " << solver_name;

  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);
  KSPDestroy(&ksp);

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0)
    << "Done cleaning up diffusion solver: " << solver_name;
}