#include"../../../ChiLua/chi_lua.h"
#include <iostream>

#include "../chi_meshhandler.h"

/** \defgroup LuaMeshHandler Mesh Handler
 * \ingroup LuaMesh
*/
//#############################################################################
/** Gets an index to a surface from a SurfaceMeshCollection.
 *
\param CollectionHandle int Handle to the SurfaceMeshCollection
\param SurfaceIndex int Index of the SurfaceMesh in the collection.

\return Handle int Handle to the extracted SurfaceMesh.
\ingroup LuaMeshHandler
\author Jan*/
int chiMeshHandlerGetSurfaceFromCollection(lua_State *L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  int coll_handle = lua_tonumber(L,1);
  int surf_index  = lua_tonumber(L,2);

  try{
    chi_mesh::MeshHandler::SurfaceMeshCollection* curItem =
      cur_hndlr.surface_mesh_collections.at(coll_handle);

    chi_mesh::SurfaceMesh* surf = curItem->at(surf_index);

    cur_hndlr.surface_mesh_stack.push_back(surf);
    int index = (int)cur_hndlr.surface_mesh_stack.size()-1;
    lua_pushnumber(L,index);

    std::cout << "chiMeshHandlerGetSurfaceFromCollection: Surface extracted = ";
    std::cout << index;
    std::cout << std::endl;
  }

  catch(const std::out_of_range& o){
    std::cerr << "ERROR: Invalid handle or index specified in "
                 "chiMeshHandlerGetSurfaceFromCollection.\n";
    exit(EXIT_FAILURE);
  }



  return 1;
}
