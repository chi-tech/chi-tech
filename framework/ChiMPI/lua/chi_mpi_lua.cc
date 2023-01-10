#include "chi_mpi_lua.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)

void chi_mpi_utils::lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiMPIBarrier);
}