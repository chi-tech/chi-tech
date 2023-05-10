#ifndef CHITECH_LBS_LUA_UTILS_H
#define CHITECH_LBS_LUA_UTILS_H

#include "chi_lua.h"

namespace lbs::common_lua_utils
{
int chiLBSSetOptions(lua_State* L);
int chiLBSSetPhiFromFieldFunction(lua_State* L);
void RegisterLuaEntities(lua_State* L);
}

#endif //CHITECH_LBS_LUA_UTILS_H
