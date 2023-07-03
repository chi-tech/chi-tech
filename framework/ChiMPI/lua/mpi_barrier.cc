#include "ChiLua/chi_lua.h"

#include "chi_runtime.h"

//#############################################################################
/** Blocks until all processes in the communicator have reached this routine.

\ingroup chiMPI
\author Jan*/
int chiMPIBarrier(lua_State *L)
{

  MPI_Barrier(Chi::mpi.comm);
  return 0;
}