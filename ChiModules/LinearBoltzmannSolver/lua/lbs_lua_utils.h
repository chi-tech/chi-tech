#ifndef LBS_LUA_UTILS_H
#define LBS_LUA_UTILS_H

#include "../lbs_linear_boltzmann_solver.h"

#include "chi_lua.h"

int chiLBSCreateSolver(lua_State *L);
int chiLBSSetProperty(lua_State *L);
int chiLBSInitialize(lua_State *L);
int chiLBSExecute(lua_State *L);
int chiLBSGetFieldFunctionList(lua_State *L);
int chiLBSGetScalarFieldFunctionList(lua_State *L);
int chiLBSWriteGroupsetAngularFlux(lua_State *L);
int chiLBSReadGroupsetAngularFlux(lua_State *L);
int chiLBSWriteFluxMoments(lua_State *L);
int chiLBSCreateAndWriteSourceMoments(lua_State *L);
int chiLBSReadFluxMomentsAndMakeSourceMoments(lua_State *L);
int chiLBSReadSourceMoments(lua_State *L);
int chiLBSReadFluxMoments(lua_State *L);
int chiLBSComputeBalance(lua_State *L);

int chiLBSCreateGroupset(lua_State *L);
int chiLBSCreateGroup(lua_State *L);
int chiLBSGroupsetAddGroups(lua_State *L);
int chiLBSGroupsetSetQuadrature(lua_State *L);
int chiLBSGroupsetSetAngleAggregationType(lua_State *L);
int chiLBSGroupsetSetAngleAggDiv(lua_State *L);
int chiLBSGroupsetSetGroupSubsets(lua_State *L);
int chiLBSGroupsetSetIterativeMethod(lua_State *L);
int chiLBSGroupsetSetResidualTolerance(lua_State *L);
int chiLBSGroupsetSetMaxIterations(lua_State *L);
int chiLBSGroupsetSetGMRESRestartIntvl(lua_State *L);
int chiLBSGroupsetSetEnableSweepLog(lua_State *L);
int chiLBSGroupsetSetWGDSA(lua_State *L);
int chiLBSGroupsetSetTGDSA(lua_State *L);

int chiLBSAddPointSource(lua_State *L);

namespace lbs
{
  namespace lua_utils
  {
    //###################################################################
    /** Obtains a pointer to a LinearBoltzmann::Solver object or an object
     * derived from LinearBoltzmann::Solver
     *
     * \param handle int Index in the chi_physics_handler where the solve object
     *                   should be located.
     * \param calling_function_name string The string used to print error messages,
     *                              should uniquely identify the calling function.
     *
     */
    lbs::SteadySolver* GetSolverByHandle(int handle, const std::string& calling_function_name);

    void RegisterLuaEntities(lua_State* L);
  }
}

#endif
