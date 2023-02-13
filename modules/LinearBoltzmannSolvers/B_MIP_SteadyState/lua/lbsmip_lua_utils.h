#ifndef CHITECH_LBSMIP_LUA_UTILS_H
#define CHITECH_LBSMIP_LUA_UTILS_H

#include "chi_lua.h"

namespace lbs::mip_steady_state_lua_utils
{
  int chiLBSMIPCreateSolver(lua_State *L);
  void RegisterLuaEntities(lua_State* L);
}//namespace lbs::lbsmip_lua_utils

#endif //CHITECH_LBSMIP_LUA_UTILS_H
