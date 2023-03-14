#include "chi_modules_lua.h"

#include "DiffusionSolver/lua/diffusion_lua.h"
#include "CFEMDiffusion/lua/ds_lua_utils.h"
#include "DFEMDiffusion/lua/ip_lua_utils.h"
#include "FVDiffusion/lua/ds_lua_utils.h"
#include "MGDiffusion/lua/mgds_lua_utils.h"
#include "A_LBSSolver/lua/lbs_lua_utils.h"
#include "Ca_DO_SteadyState/lua/lbs_disc_ord_ss_lua_utils.h"
#include "Cc_DO_KEigenvalue/lua/lbkes_lua_utils.h"
#include "Cd_DO_SteadyState_Adjoint/lua/lbsadj_lua_utils.h"
#include "D_DO_RZ_SteadyState/lua/lbs_curvilinear_solver_lua.h"
#include "D_DO_Transient/lua/lbts_lua_utils.h"
#include "B_MIP_SteadyState/lua/lbsmip_lua_utils.h"

void chi_modules::lua_utils::RegisterLuaEntities(lua_State *L)
{
  lbs::common_lua_utils::RegisterLuaEntities(L);
  lbs::disc_ord_steady_state_lua_utils::RegisterLuaEntities(L);
  diffusion_solver::lua_utils::RegisterLuaEntities(L);
  lbs_curvilinear::lua_utils::RegisterLuaEntities(L);
  lbs::k_eigenvalue_lua_utils::RegisterLuaEntities(L);
  lbs::adjoint_lua_utils::RegisterLuaEntities(L);
  cfem_diffusion::cfem_diffusion_lua_utils::RegisterLuaEntities(L);
  dfem_diffusion::dfem_diffusion_lua_utils::RegisterLuaEntities(L);
  lbs::lbts_lua_utils::RegisterLuaEntities(L);
  mg_diffusion::mgd_lua_utils::RegisterLuaEntities(L);
  fv_diffusion::fv_diffusion_lua_utils::RegisterLuaEntities(L);
  lbs::mip_steady_state_lua_utils::RegisterLuaEntities(L);
}