#ifndef CFEM_DIFFUSION_LUA_UTILS_H
#define CFEM_DIFFUSION_LUA_UTILS_H

#include"ChiLua/chi_lua.h"
#include "../cfem_diffusion_solver.h"

int chiCFEMDiffusionSolverCreate(lua_State *L);
int chiCFEMDiffusionSetBCProperty(lua_State *L);


namespace cfem_diffusion
{
  namespace cfem_diffusion_lua_utils
  {
    //###################################################################
    /** Obtains a pointer to a cfem_diffusion::Solver object or an object
     * derived from chi_physics::Solver
     *
     * \param handle int Index in the chi_physics_handler where the solve object
     *                   should be located.
     * \param calling_function_name string The string used to print error messages,
     *                              should uniquely identify the calling function.
     *
     */
    cfem_diffusion::Solver& GetSolverByHandle(
      int handle, const std::string& calling_function_name);

    void RegisterLuaEntities(lua_State *L);
  }//namespace cfem_diffusion_lua_utils
}//namespace cfem_diffusion


#endif //CFEM_DIFFUSION_LUA_UTILS_H