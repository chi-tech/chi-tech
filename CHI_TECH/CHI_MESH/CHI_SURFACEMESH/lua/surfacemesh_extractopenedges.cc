#include"../../../CHI_LUA/chi_lua.h"

#include <iostream>
#include <sstream>
#include "../chi_surfacemesh.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"

#include <chi_log.h>

extern CHI_LOG chi_log;


//#############################################################################
/** Exports all open edges of a surface mesh to file. This is used mostly
 * for graphical error checking.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.
\param FileName char Filename to which the edges are to be exported.

\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshExtractOpenEdgesToObj(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("chiSurfaceMeshExtractOpenEdgesToObj",2,num_args);

  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int         surf_handle = lua_tonumber(L,1);
  const char* file_name   = lua_tostring(L,2);

  try{
    chi_mesh::SurfaceMesh* curItem =
      cur_hndlr->surface_mesh_stack.at(surf_handle);

    curItem->ExtractOpenEdgesToObj(file_name);
  }

  catch(std::out_of_range o){
    std::cerr << "ERROR: Invalid index to surface mesh.\n";
    exit(EXIT_FAILURE);
  }
  return 0;
}