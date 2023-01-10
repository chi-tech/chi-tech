#include "lbsadj_lua_utils.h"

void lbs_adjoint::lua_utils::RegisterLuaEntities(lua_State* L)
{
  lua_register(L, "chiAdjointSolverCreate",
               lbs_adjoint::lua_utils::chiAdjointSolverCreate);
  lua_register(L, "chiAdjointSolverAddResponseFunction",
               lbs_adjoint::lua_utils::chiAdjointSolverAddResponseFunction);
  lua_register(L, "chiAdjointSolverMakeExpRepFromP1Moments",
               lbs_adjoint::lua_utils::chiAdjointSolverMakeExpRepFromP1Moments);
  lua_register(L, "chiAdjointSolverExportImportanceMapBinary",
               lbs_adjoint::lua_utils::chiAdjointSolverExportImportanceMapBinary);
  lua_register(L, "chiAdjointSolverComputeInnerProduct",
               lbs_adjoint::lua_utils::chiAdjointSolverComputeInnerProduct);
  lua_register(L, "chiAdjointSolverReadFluxMomentsToBuffer",
               lbs_adjoint::lua_utils::chiAdjointSolverReadFluxMomentsToBuffer);
  lua_register(L, "chiAdjointSolverApplyFluxMomentBuffer",
               lbs_adjoint::lua_utils::chiAdjointSolverApplyFluxMomentBuffer);
}