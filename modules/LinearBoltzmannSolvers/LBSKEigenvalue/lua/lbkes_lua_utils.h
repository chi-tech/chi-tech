#ifndef LBKES_LUA_UTILS_H
#define LBKES_LUA_UTILS_H

#include "../lbkes_k_eigenvalue_solver.h"

#include "chi_lua.h"

int chiLBKESCreateSolver(lua_State *L);
int chiLBKESSetProperty(lua_State *L);

namespace lbs::k_eigenvalue_lua_utils
{
  void RegisterLuaEntities(lua_State *L);
}//namespace lbs::k_eigenvalue_lua_utils

#endif //LBKES_LUA_UTILS_H
