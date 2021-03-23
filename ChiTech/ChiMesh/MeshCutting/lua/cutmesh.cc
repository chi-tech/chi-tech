#include "ChiLua/chi_lua.h"

#include "../meshcutting.h"
#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"

//###################################################################
/**Cuts a mesh.*/
int chiCutMesh(lua_State* L)
{
  auto fname = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname,L,1);
  LuaCheckNilValue(fname,L,2);

  std::vector<double> p_raw;
  std::vector<double> n_raw;
  LuaPopulateVectorFrom1DArray(fname,L,1,p_raw);
  LuaPopulateVectorFrom1DArray(fname,L,2,n_raw);

  auto handler = chi_mesh::GetCurrentHandler();

  auto grid = handler->GetGrid();

  chi_mesh::Vector3 p(p_raw[0],p_raw[1],p_raw[2]);
  chi_mesh::Vector3 n(n_raw[0],n_raw[1],n_raw[2]);

  std::cout << p.PrintS() << std::endl;
  std::cout << n.PrintS() << std::endl;
  chi_mesh::mesh_cutting::CutMeshWithPlane(*grid,p,n);

  return 0;
}