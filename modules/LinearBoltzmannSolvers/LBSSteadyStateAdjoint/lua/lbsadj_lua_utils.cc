#include "lbsadj_lua_utils.h"

void lbs::lua_utils::RegisterLuaEntities(lua_State* L)
{
  lua_register(L, "chiAdjointSolverCreate",
               lbs::lua_utils::chiAdjointSolverCreate);
  lua_register(L, "chiAdjointSolverAddResponseFunction",
               lbs::lua_utils::chiAdjointSolverAddResponseFunction);
  lua_register(L, "chiAdjointSolverMakeExpRepFromP1Moments",
               lbs::lua_utils::chiAdjointSolverMakeExpRepFromP1Moments);
  lua_register(L, "chiAdjointSolverExportImportanceMapBinary",
               lbs::lua_utils::chiAdjointSolverExportImportanceMapBinary);
  lua_register(L, "chiAdjointSolverComputeInnerProduct",
               lbs::lua_utils::chiAdjointSolverComputeInnerProduct);
  lua_register(L, "chiAdjointSolverReadFluxMomentsToBuffer",
               lbs::lua_utils::chiAdjointSolverReadFluxMomentsToBuffer);
  lua_register(L, "chiAdjointSolverApplyFluxMomentBuffer",
               lbs::lua_utils::chiAdjointSolverApplyFluxMomentBuffer);
}