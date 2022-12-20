#ifndef LBSADJOINTSOLVER_LUA_UTILS_H
#define LBSADJOINTSOLVER_LUA_UTILS_H

#include "LBSAdjointSolver/lbsadj_solver.h"

#include "chi_lua.h"

namespace lbs_adjoint::lua_utils
{
  int chiAdjointSolverCreate(lua_State* L);
  int chiAdjointSolverAddResponseFunction(lua_State* L);
  int chiAdjointSolverMakeExpRepFromP1Moments(lua_State* L);
  int chiAdjointSolverExportImportanceMapBinary(lua_State* L);
  int chiAdjointSolverComputeInnerProduct(lua_State* L);

  int chiAdjointSolverReadFluxMomentsToBuffer(lua_State* L);
  int chiAdjointSolverApplyFluxMomentBuffer(lua_State* L);

  void RegisterLuaEntities(lua_State* L);
}//namespace lbs_adjoint

#endif //LBSADJOINTSOLVER_LUA_UTILS_H