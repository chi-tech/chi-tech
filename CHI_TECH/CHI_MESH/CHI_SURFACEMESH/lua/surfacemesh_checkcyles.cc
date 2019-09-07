#include"../../../CHI_LUA/chi_lua.h"

#include <iostream>
#include <sstream>
#include "../chi_surfacemesh.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"

#include <chi_log.h>

extern ChiLog chi_log;


//#############################################################################
/** Builds sweep ordering for a number of angles and checks whether any
 * cyclic dependencies are encountered.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.
\param NumAngles int Number of azimuthal angles to use for checking cycles.

\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshCheckCycles(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("chiSurfaceMeshCheckCycles",2,num_args);

  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int surf_handle = lua_tonumber(L,1);
  int num_angles  = lua_tonumber(L,2);

  try{
    chi_mesh::SurfaceMesh* curItem =
      cur_hndlr->surface_mesh_stack.at(surf_handle);

    curItem->CheckCyclicDependencies(num_angles);
  }

  catch(std::out_of_range o){
    std::cerr << "ERROR: Invalid index to surface mesh.\n";
    exit(EXIT_FAILURE);
  }
  return 0;
}