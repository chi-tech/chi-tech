#include "../../ChiLua/chi_lua.h"
#include <iostream>

#include "../chi_mpi.h"

extern ChiMPI& chi_mpi;

/** \defgroup chiMPI E MPI Utilities

## Lua available variables

- *chi_location_id* - (int) Process number for current process
- *chi_number_of_processes* - (int) Total number of processes
 * */

//#############################################################################
/** Blocks until all processes in the communicator have reached this routine.

\ingroup chiMPI
\author Jan*/
int chiMPIBarrier(lua_State *L)
{

  MPI_Barrier(MPI_COMM_WORLD);
  return 0;
}