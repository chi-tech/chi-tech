#include "k_eigenvalue_solver.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

#include <iomanip>


using namespace LinearBoltzmann;

//###################################################################
void KEigenvalue::Solver::Execute()
{
  MPI_Barrier(MPI_COMM_WORLD);

  PowerIteration();
  chi_log.Log(LOG_0) << "KEigenvalue::Solver execution completed\n";

  MPI_Barrier(MPI_COMM_WORLD);
}