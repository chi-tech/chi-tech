#include"../../../ChiLua/chi_lua.h"

#include <iostream>
#include <sstream>
#include "../chi_surfacemesh.h"
#include "../../MeshHandler/chi_meshhandler.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//#############################################################################
/** Splits a SurfaceMesh by patch.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.

\return Handle int Handle to the patch collection.
\return Count int Number of patches found.
\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshSplitByPatch(lua_State *L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  int surf_handle = lua_tonumber(L,1);

  try{
    chi_mesh::SurfaceMesh* curItem =
      cur_hndlr.surface_mesh_stack.at(surf_handle);

    chi_mesh::MeshHandler::SurfaceMeshCollection* new_coll =
      new chi_mesh::MeshHandler::SurfaceMeshCollection;

    curItem->SplitByPatch(*new_coll);

    cur_hndlr.surface_mesh_collections.push_back(new_coll);
    int index = cur_hndlr.surface_mesh_collections.size()-1;
    lua_pushnumber(L,index);
    lua_pushnumber(L,new_coll->size());

    std::stringstream outtext;
    outtext << "chiSurfaceMeshSplitByPatch: Number of patches found = ";
    outtext << new_coll->size();
    outtext << std::endl;

    chi_log.Log(LOG_ALLVERBOSE_2) << outtext.str();
  }

  catch(const std::out_of_range& o){
    std::cerr << "ERROR: Invalid index to surface mesh.\n";
    exit(EXIT_FAILURE);
  }




  return 2;
}

