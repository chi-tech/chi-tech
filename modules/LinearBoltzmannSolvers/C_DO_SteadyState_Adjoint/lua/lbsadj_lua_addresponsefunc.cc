#include "C_DO_SteadyState_Adjoint/lbsadj_solver.h"
#include "lbsadj_lua_utils.h"

#include "chi_runtime.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"


namespace lbs::adjoint_lua_utils
{

int chiAdjointSolverAddResponseFunction(lua_State* L)
{
  const std::string fname = "chiAdjointSolverAddResponseFunction";
  const int num_args = lua_gettop(L);
  if (num_args < 3)
    LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  LuaCheckIntegerValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);
  LuaCheckNumberValue(fname, L, 3);

  const int solver_handle     = lua_tointeger(L, 1);
  const std::string qoi_name = lua_tostring(L,2);
  const int logvol_handle    = lua_tointeger(L,3);

  std::string lua_function;
  if (num_args == 4)
  {
    LuaCheckNilValue(fname, L, 4);

    lua_function = lua_tostring(L,4);
  }

  auto& solver = chi::GetStackItem<lbs::DiscOrdSteadyStateAdjointSolver>(
    chi::solver_stack, solver_handle, fname);

  auto p_logical_volume = chi::GetStackItemPtr(
    chi::logicvolume_stack, logvol_handle, fname);

  size_t qoi_index = solver.AddResponseFunction(qoi_name,
                                                 p_logical_volume,
                                                 lua_function);
  lua_pushinteger(L, static_cast<lua_Integer>(qoi_index));

  return 1;
}

}//namespace lbs

