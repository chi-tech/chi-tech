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
    void RegisterLuaEntities(lua_State *L);
  }//namespace cfem_diffusion_lua_utils
}//namespace cfem_diffusion


#endif //CFEM_DIFFUSION_LUA_UTILS_H