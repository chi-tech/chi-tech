#include "ChiLua/chi_lua.h"

#include "../meshcutting.h"
#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"

//###################################################################
/**Cuts a mesh.*/
int chiCutMesh(lua_State* L)
{
  auto handler = chi_mesh::GetCurrentHandler();

  auto grid = handler->GetGrid();

  chi_mesh::Vector3 p(0.0,0.0,0.0);
  chi_mesh::Vector3 n(1.0,0.0,0.0);
  chi_mesh::mesh_cutting::CutMeshWithPlane(*grid,p,n);

  return 0;
}