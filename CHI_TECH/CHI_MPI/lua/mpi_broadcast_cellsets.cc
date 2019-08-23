#include "../../CHI_LUA/chi_lua.h"
#include <iostream>

#include "../chi_mpi.h"

extern CHI_MPI chi_mpi;

/** \defgroup chiMPI E MPI Utilities

## Lua available variables

- *chi_location_id* - (int) Process number for current process
- *chi_number_of_processes* - (int) Total number of processes
 * */

//#############################################################################
/** Broadcasts mesh to all child processes.

\ingroup chiMPI
\author Jan*/
int chiMPIBroadcastCellsets(lua_State *L)
{

  chi_mpi.BroadcastCellSets();
  return 0;
}