#include"../../../CHI_LUA/chi_lua.h"

#include <iostream>
#include <sstream>
#include "../chi_surfacemesh.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include <chi_log.h>

extern CHI_LOG chi_log;



//############################################################################# Create
/** Exports mesh as a .poly format.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.
\param FileName char* Path to the file to be imported.

\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshExportPolyFile(lua_State* L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("chiSurfaceMeshExportPolyFile",2, num_args);

  int handle = lua_tonumber(L,1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);



  try{
    chi_mesh::SurfaceMesh* curItem = cur_hndlr->surface_mesh_stack.at(handle);

    curItem->ExportToPolyFile(temp);
  }

  catch(std::out_of_range o){
    std::cerr << "ERROR: Invalid index to surface mesh.\n";
    exit(EXIT_FAILURE);
  }
  return 0;
}