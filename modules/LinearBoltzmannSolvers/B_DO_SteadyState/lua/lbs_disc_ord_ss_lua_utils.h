#ifndef LBS_LUA_UTILS_H
#define LBS_LUA_UTILS_H

#include "chi_lua.h"

namespace lbs::disc_ord_steady_state_lua_utils
{
  int chiLBSCreateSolver(lua_State* L);
  int chiLBSComputeBalance(lua_State* L);
  int chiLBSComputeLeakage(lua_State* L);

  void RegisterLuaEntities(lua_State* L);
}

#endif
