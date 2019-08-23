#include "../../CHI_LUA/chi_lua.h"
#include <iostream>

#include "../chi_mpi.h"

extern CHI_MPI chi_mpi;

//#############################################################################
/** Receive mesh from parent process.

\ingroup chiMPI
\author Jan*/
int chiMPIReceiveCellsets(lua_State *L)
{
  chi_mpi.ReceiveCellSets();
  return 0;
}