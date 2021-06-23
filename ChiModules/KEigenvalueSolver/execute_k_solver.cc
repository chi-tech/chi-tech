#include "k_eigenvalue_solver.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI&      chi_mpi;
extern ChiLog&     chi_log;

#include <iomanip>


using namespace LinearBoltzmann;

//###################################################################
void KEigenvalue::Solver::ExecuteKSolver()
{
  MPI_Barrier(MPI_COMM_WORLD);

  for (auto& groupset : group_sets) {
    ComputeSweepOrderings(groupset);
    InitFluxDataStructures(groupset);
  }

  PowerIteration();

  for (auto& groupset : group_sets) {
    ResetSweepOrderings(groupset);
  }


  chi_log.Log(LOG_0) << "KEigenvalueSolver execution completed\n";

  MPI_Barrier(MPI_COMM_WORLD);
}