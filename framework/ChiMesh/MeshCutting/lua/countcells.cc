#include "ChiLua/chi_lua.h"

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMesh/LogicalVolume/LogicalVolume.h"

#include "chi_runtime.h"

#include "chi_log.h"

//###################################################################
/**Counts the number of cells with a logical volume.*/
int chiCountMeshInLogicalVolume(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  //======================================== Arg checking
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(__FUNCTION__,1,num_args);

  LuaCheckNilValue(__FUNCTION__,L,1);

  int log_vol_handle = lua_tonumber(L,1);

  auto& handler = chi_mesh::GetCurrentHandler();

  const auto& log_vol = Chi::GetStackItem<chi_mesh::LogicalVolume>(
    Chi::object_stack, log_vol_handle, fname);

  auto& grid = handler.GetGrid();

  size_t count = grid->CountCellsInLogicalVolume(log_vol);

  lua_pushinteger(L,int(count));
  return 1;
}