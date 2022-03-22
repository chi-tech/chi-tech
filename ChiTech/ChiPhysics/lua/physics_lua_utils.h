#ifndef CHI_PHYSICS_LUA_UTILS_H
#define CHI_PHYSICS_LUA_UTILS_H

#include "ChiPhysics/SolverBase/chi_solver.h"
#include "chi_lua.h"

int chiSolverAddRegion(lua_State *L);
int chiSolverInitialize(lua_State *L);
int chiSolverExecute(lua_State *L);
int chiSolverSetBasicOption(lua_State *L);

int chiPhysicsAddMaterial(lua_State *L);
int chiPhysicsMaterialAddProperty(lua_State *L);
int chiPhysicsMaterialSetProperty(lua_State *L);
int chiPhysicsMaterialGetProperty(lua_State *L);
int chiPhysicsMaterialModifyTotalCrossSection(lua_State *L);


namespace chi_physics
{
  namespace lua_utils
  {
    //###################################################################
    /** Obtains a pointer to a chi_physics::Solver object or an object
     * derived from chi_physics::Solver
     *
     * \param handle int Index in the chi_physics_handler where the solve object
     *                   should be located.
     * \param calling_function_name string The string used to print error messages,
     *                              should uniquely identify the calling function.
     *
     */
     chi_physics::Solver* GetSolverByHandle(int handle, const std::string& calling_function_name);

     void RegisterLuaEntities(lua_State* L);
  }
}

#endif //CHI_PHYSICS_LUA_UTILS_H