#ifndef CFEM_DIFFUSION_LUA_UTILS_H
#define CFEM_DIFFUSION_LUA_UTILS_H

#include"ChiLua/chi_lua.h"

int chiCFEMDiffusionSolverCreate(lua_State *L);
int chiCFEMDiffusionSetBCProperty(lua_State *L);


#endif //CFEM_DIFFUSION_LUA_UTILS_H