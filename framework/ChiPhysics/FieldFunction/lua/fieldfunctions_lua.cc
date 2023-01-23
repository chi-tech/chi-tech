#include "fieldfunctions_lua.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)

void chi_physics::field_function_lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiGetFieldFunctionHandleByName);
  LUA_FMACRO1(chiExportFieldFunctionToVTK);
  LUA_FMACRO1(chiExportMultiFieldFunctionToVTK);
}