#include"../../../ChiLua/chi_lua.h"

#include <iostream>
#include <sstream>
#include "../chi_surfacemesh.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"

#include <chi_log.h>

extern ChiLog chi_log;

//#############################################################################
/** Gets a list of edge loops for the given surface mesh.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.

\return Handle int Handle to the edge loops.
\return Count int Number of edge loops found.
\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshGetEdgeLoops(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int surf_handle = lua_tonumber(L,1);

  try{
    chi_mesh::SurfaceMesh* curItem =
            cur_hndlr->surface_mesh_stack.at(surf_handle);
    chi_mesh::EdgeLoopCollection* loops;
    loops = curItem->GetEdgeLoops();

    cur_hndlr->edge_loop_collections.push_back(loops);
    int index = cur_hndlr->edge_loop_collections.size()-1;
    lua_pushnumber(L,index);
    lua_pushnumber(L,loops->size());

    std::stringstream outtext;
    outtext << "chiSurfaceMeshGetEdgeLoops: Number of open edge loops found = ";
    outtext << loops->size();
    outtext << std::endl;

    chi_log.Log(LOG_ALLVERBOSE_2) << outtext.str();
  }

  catch(std::out_of_range o){
    std::cerr << "ERROR: Invalid index to surface mesh.\n";
    exit(EXIT_FAILURE);
  }




  return 2;
}

//#############################################################################
/** Gets a list of edge loops for the given surface mesh's polygon faces.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.

\return Handle int Handle to the edge loops.
\return Count int Number of edge loops found.
\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshGetEdgeLoopsPoly(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int surf_handle = lua_tonumber(L,1);

  try{
    chi_mesh::SurfaceMesh* curItem =
      cur_hndlr->surface_mesh_stack.at(surf_handle);
    chi_mesh::EdgeLoopCollection* loops;
    loops = curItem->GetEdgeLoopsPoly();

    cur_hndlr->edge_loop_collections.push_back(loops);
    int index = cur_hndlr->edge_loop_collections.size()-1;
    lua_pushnumber(L,index);
    lua_pushnumber(L,loops->size());

    std::stringstream outtext;
    outtext << "chiSurfaceMeshGetEdgeLoops: Number of open edge loops found = ";
    outtext << loops->size();
    outtext << std::endl;

    chi_log.Log(LOG_ALLVERBOSE_2) << outtext.str();
  }

  catch(std::out_of_range o){
    std::cerr << "ERROR: Invalid index to surface mesh.\n";
    exit(EXIT_FAILURE);
  }




  return 2;
}