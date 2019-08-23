#include "../../CHI_LUA/chi_lua.h"
#include <iostream>

#include "../chi_mpi.h"

extern CHI_MPI chi_mpi;

//#############################################################################
/** Blocks until all processes in the communicator have reached this routine.

\ingroup chiMPI
\author Jan*/
int chiMPIBarrier(lua_State *L)
{

  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}