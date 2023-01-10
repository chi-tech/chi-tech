#include "lbs_linear_boltzmann_solver.h"

#include <chi_mpi.h>

//###################################################################
/** Initialize the solver.*/
void lbs::SteadySolver::Initialize()
{
  PerformInputChecks(); //a
  PrintSimHeader(); //b
  MPI_Barrier(MPI_COMM_WORLD);

  //================================================== Initialize materials
  InitMaterials(); //c

  //================================================== Init spatial discretization
  InitializeSpatialDiscretization(); //d

  //================================================== Initialize groupsets
  InitializeGroupsets(); //e

  //================================================== Compute n. moments
  ComputeNumberOfMoments(); //f

  //================================================== Initialize parrays
  InitializeParrays();//g

  //================================================== Initialize boundaries
  InitializeBoundaries();

  //================================================== Initialize sources
  InitializePointSources();

  //================================================== Initialize source func
  using namespace std::placeholders;
  active_set_source_function =
    std::bind(&SteadySolver::SetSource, this, _1, _2, _3);

}
