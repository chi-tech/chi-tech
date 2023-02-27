#include "lbsadj_lua_utils.h"

void lbs::adjoint_lua_utils::RegisterLuaEntities(lua_State* L)
{
  lua_register(L, "chiAdjointSolverCreate",
               lbs::adjoint_lua_utils::chiAdjointSolverCreate);
  lua_register(L, "chiAdjointSolverAddResponseFunction",
               lbs::adjoint_lua_utils::chiAdjointSolverAddResponseFunction);
  lua_register(L, "chiAdjointSolverMakeExpRepFromP1Moments",
               lbs::adjoint_lua_utils::chiAdjointSolverMakeExpRepFromP1Moments);
  lua_register(L, "chiAdjointSolverExportImportanceMapBinary",
               lbs::adjoint_lua_utils::chiAdjointSolverExportImportanceMapBinary);
  lua_register(L, "chiAdjointSolverComputeInnerProduct",
               lbs::adjoint_lua_utils::chiAdjointSolverComputeInnerProduct);
  lua_register(L, "chiAdjointSolverReadFluxMomentsToBuffer",
               lbs::adjoint_lua_utils::chiAdjointSolverReadFluxMomentsToBuffer);
  lua_register(L, "chiAdjointSolverApplyFluxMomentBuffer",
               lbs::adjoint_lua_utils::chiAdjointSolverApplyFluxMomentBuffer);
}