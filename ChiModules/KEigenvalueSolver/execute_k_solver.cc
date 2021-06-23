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
  int gs = 0;

  chi_log.Log(LOG_0)
    << "\n********* Initializing Groupset " << gs << std::endl;

  group_sets[gs].BuildDiscMomOperator(options.scattering_order,
                                      options.geometry_type);
  group_sets[gs].BuildMomDiscOperator(options.scattering_order,
                                      options.geometry_type);
  group_sets[gs].BuildSubsets();

  ComputeSweepOrderings(group_sets[gs]);
  InitFluxDataStructures(group_sets[gs]);

  PowerIteration();

  ResetSweepOrderings(group_sets[gs]);
 
  chi_log.Log(LOG_0) << "KEigenvalueSolver execution completed\n";

  MPI_Barrier(MPI_COMM_WORLD);
}