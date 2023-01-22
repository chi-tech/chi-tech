#ifndef CHITECH_FIELDFUNCTIONS_LUA_H
#define CHITECH_FIELDFUNCTIONS_LUA_H

#include "chi_lua.h"

int chiGetFieldFunctionHandleByName(lua_State *L);
int chiExportFieldFunctionToVTK(lua_State *L);
int chiExportMultiFieldFunctionToVTK(lua_State *L);

namespace chi_physics::field_function_lua_utils
{
  void RegisterLuaEntities(lua_State *L);
}//namespace chi_physics::field_function_lua_utils


#endif //CHITECH_FIELDFUNCTIONS_LUA_H
