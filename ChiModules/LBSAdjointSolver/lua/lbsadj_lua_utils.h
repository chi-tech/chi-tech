#ifndef LBSADJOINTSOLVER_LUA_UTILS_H
#define LBSADJOINTSOLVER_LUA_UTILS_H

#include "LBSAdjointSolver/lbsadj_solver.h"

#include "chi_lua.h"

namespace lbs_adjoint::lua_utils
  {
    //###################################################################
    /** Obtains a pointer to a lbs_adjoint::AdjointSolver object or an object
     * derived from lbs_adjoint::AdjointSolver
     *
     * \param handle int Index in the chi_physics_handler where the solve object
     *                   should be located.
     * \param calling_function_name string The string used to print error messages,
     *                              should uniquely identify the calling function.
     *
     */
    lbs_adjoint::AdjointSolver& GetSolverByHandle(
      int handle, const std::string& calling_function_name);

    int chiAdjointSolverCreate(lua_State* L);
    int chiAdjointSolverAddResponseFunction(lua_State* L);
    int chiAdjointSolverMakeExpRepFromP1Moments(lua_State* L);
    int chiAdjointSolverExportImportanceMapBinary(lua_State* L);
    int chiAdjointSolverComputeInnerProduct(lua_State* L);

    void RegisterLuaEntities(lua_State* L);
  }//namespace lbs_adjoint

#endif //LBSADJOINTSOLVER_LUA_UTILS_H