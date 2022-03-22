#include "chi_modules_lua.h"

#include "LinearBoltzmannSolver/lua/lbs_lua_utils.h"
#include "DiffusionSolver/lua/diffusion_lua.h"
#include "LBSCurvilinear/lua/lbs_curvilinear_solver_lua.h"
#include "LBKEigenvalueSolver/lua/lbkes_lua_utils.h"

void chi_modules::lua_utils::RegisterLuaEntities(lua_State *L)
{
  LinearBoltzmann::lua_utils::RegisterLuaEntities(L);
  diffusion_solver::lua_utils::RegisterLuaEntities(L);
  LBSCurvilinear::lua_utils::RegisterLuaEntities(L);
  LinearBoltzmann::k_eigenvalue_lua_utils::RegisterLuaEntities(L);
}