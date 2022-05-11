#include "LBSAdjointSolver/lbsadj_solver.h"
#include "lbsadj_lua_utils.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

namespace lbs_adjoint
{
namespace lua_utils
{

int chiAdjointSolverExportImportanceMapBinary(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  LuaCheckIntegerValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);

  const int solver_index     = lua_tointeger(L,1);
  const std::string file_name = lua_tostring(L,2);

  auto& solver = lbs_adjoint::lua_utils::GetSolverByHandle(solver_index,fname);

  solver.ExportImportanceMap(file_name);

  return 0;
}

}//namespace lua_utils
}//namespace lbs_adjoint

