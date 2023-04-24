#ifndef CHITECH_PHYSICS_SOLVER_LUA_UTILS_H
#define CHITECH_PHYSICS_SOLVER_LUA_UTILS_H

#include "chi_lua.h"
namespace chi_physics::lua_utils
{
int chiSolverCreate(lua_State* L);
int chiSolverInitialize(lua_State* L);
int chiSolverExecute(lua_State* L);
int chiSolverStep(lua_State* L);
int chiSolverAdvance(lua_State* L);
int chiSolverSetBasicOption(lua_State* L);
int chiSolverGetName(lua_State* L);
int chiSolverGetFieldFunctionList(lua_State* L);
} // namespace chi_physics::lua_utils

#endif // CHITECH_PHYSICS_SOLVER_LUA_UTILS_H
