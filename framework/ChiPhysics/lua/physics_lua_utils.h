#ifndef CHI_PHYSICS_LUA_UTILS_H
#define CHI_PHYSICS_LUA_UTILS_H

#include "ChiPhysics/SolverBase/chi_solver.h"
#include "chi_lua.h"

int chiSolverInitialize(lua_State* L);
int chiSolverExecute(lua_State* L);
int chiSolverStep(lua_State* L);
int chiSolverSetBasicOption(lua_State* L);
int chiSolverGetName(lua_State* L);
int chiSolverGetFieldFunctionList(lua_State* L);

int chiPhysicsAddMaterial(lua_State* L);
int chiPhysicsMaterialAddProperty(lua_State* L);
int chiPhysicsMaterialSetProperty(lua_State* L);
int chiPhysicsMaterialGetProperty(lua_State* L);
int chiPhysicsMaterialModifyTotalCrossSection(lua_State* L);


namespace chi_physics::lua_utils
{
   void RegisterLuaEntities(lua_State* L);
}

#endif //CHI_PHYSICS_LUA_UTILS_H