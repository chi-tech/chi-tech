#include "ChiLua/chi_lua.h"

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
extern ChiLog& chi_log;

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

  chi_mesh::LogicalVolume* log_vol;
  try {log_vol = handler.logicvolume_stack.at(log_vol_handle);}
  catch (const std::out_of_range& oor)
  {throw std::invalid_argument(fname + ": Invalid handle to logical volume.");}

  auto grid = handler.GetGrid();

  size_t count = grid->CountCellsInLogicalVolume(*log_vol);

  lua_pushinteger(L,int(count));
  return 1;
}