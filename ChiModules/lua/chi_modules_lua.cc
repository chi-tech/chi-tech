#include "chi_modules_lua.h"

#include "LinearBoltzmannSolver/lua/lbs_lua_utils.h"
#include "DiffusionSolver/lua/diffusion_lua.h"
#include "LBSCurvilinear/lua/lbs_curvilinear_solver_lua.h"
#include "LBKEigenvalueSolver/lua/lbkes_lua_utils.h"
#include "LBSAdjointSolver/lua/lbsadj_lua_utils.h"
#include "CFEMDiffusion/lua/ds_lua_utils.h"
#include "DFEMDiffusion/lua/ip_lua_utils.h"
#include "MGDiffusion/lua/mgds_lua_utils.h"

void chi_modules::lua_utils::RegisterLuaEntities(lua_State *L)
{
  lbs::lua_utils::RegisterLuaEntities(L);
  diffusion_solver::lua_utils::RegisterLuaEntities(L);
  lbs_curvilinear::lua_utils::RegisterLuaEntities(L);
  lbs::k_eigenvalue_lua_utils::RegisterLuaEntities(L);
  lbs_adjoint::lua_utils::RegisterLuaEntities(L);
  cfem_diffusion::cfem_diffusion_lua_utils::RegisterLuaEntities(L);
  dfem_diffusion::dfem_diffusion_lua_utils::RegisterLuaEntities(L);
  dfem_diffusion::dfem_diffusion_lua_utils::RegisterLuaEntities(L);
  mg_diffusion::mgd_lua_utils::RegisterLuaEntities(L);
}