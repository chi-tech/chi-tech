#include "LBSAdjointSolver/lbsadj_solver.h"
#include "lbsadj_lua_utils.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

namespace lbs_adjoint
{
namespace lua_utils
{

int chiAdjointSolverSetQOI(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 3)
    LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  LuaCheckIntegerValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);
  LuaCheckNumberValue(fname, L, 3);

  const int solver_index     = lua_tointeger(L,1);
  const std::string qoi_name = lua_tostring(L,2);
  const int logvol_handle    = lua_tointeger(L,3);

  auto solver = lbs_adjoint::lua_utils::GetSolverByHandle(solver_index,fname);

  auto mesh_handler = chi_mesh::GetCurrentHandler();

  chi_mesh::LogicalVolume* logical_volume;
  try {logical_volume = mesh_handler->logicvolume_stack.at(logvol_handle);}
  catch (const std::out_of_range& oor)
  {throw std::invalid_argument(fname + ": Invalid handle to logical volume.");}

  size_t qoi_index = solver->SetQOI(qoi_name, *logical_volume, "");
  lua_pushinteger(L, static_cast<lua_Integer>(qoi_index));

  return 1;
}

}//namespace lua_utils
}//namespace lbs_adjoint

