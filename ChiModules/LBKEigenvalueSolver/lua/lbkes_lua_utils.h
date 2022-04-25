#ifndef LBKES_LUA_UTILS_H
#define LBKES_LUA_UTILS_H

#include "../lbkes_k_eigenvalue_solver.h"

#include "chi_lua.h"

int chiLBKESCreateSolver(lua_State *L);
int chiLBKESInitialize(lua_State *L);
int chiLBKESExecute(lua_State *L);
int chiLBKESSetProperty(lua_State *L);

namespace lbs
{
  namespace k_eigenvalue_lua_utils
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
    lbs::KEigenvalueSolver* GetSolverByHandle(
        int handle, const std::string& calling_function_name);

    void RegisterLuaEntities(lua_State *L);
  }//namespace k_eigenvalue_lua_utils
}//namespace LinearBoltzmann

#endif //LBKES_LUA_UTILS_H
