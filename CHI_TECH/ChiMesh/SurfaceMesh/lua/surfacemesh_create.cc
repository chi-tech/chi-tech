#include"../../../ChiLua/chi_lua.h"
#include<iostream>
#include "../chi_surfacemesh.h"
#include "../../MeshHandler/chi_meshhandler.h"

#include <chi_log.h>

extern ChiLog chi_log;

/** \defgroup LuaSurfaceMesh Surface Meshes
 * \ingroup LuaMesh
*/

//#############################################################################
/** Creates a new empty surface mesh.

\return Handle int Handle to the created surface mesh.
\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshCreate(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();
  chi_mesh::SurfaceMesh* new_mesh = new chi_mesh::SurfaceMesh;

  cur_hndlr->surface_mesh_stack.push_back(new_mesh);

  int index = cur_hndlr->surface_mesh_stack.size()-1;
  lua_pushnumber(L,index);

  chi_log.Log(LOG_ALLVERBOSE_2) << "chiSurfaceMeshCreate: "
                                         "Empty SurfaceMesh object, "
                                      << index << ", created" << std::endl;

  return 1;
}