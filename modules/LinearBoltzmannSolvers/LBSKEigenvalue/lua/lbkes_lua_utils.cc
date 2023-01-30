#include "lbkes_lua_utils.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)
#define LUA_CMACRO1(x,y) \
        lua_pushnumber(L, y); \
        lua_setglobal(L, #x)

void lbs::k_eigenvalue_lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiLBKESCreateSolver);
  LUA_FMACRO1(chiLBKESSetProperty);
}