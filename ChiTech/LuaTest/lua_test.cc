#include <ChiLua/chi_lua.h>

#include "lua_test.h"

#include "unit_tests.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_runtime.h"
#include "chi_mpi.h"

;

#define LUA_FMACRO1(x) lua_register(L, #x, x)

//###################################################################
/**This is a lua test function.
\param argument1 Any Argument of any type.
\ingroup LuaGeneralUtilities
 */
int chiLuaTest(lua_State* L)
{
  chi::log.Log() << "Hello from chiLuaTest(). Process count: "
                << chi::mpi.process_count;
  chi::log.LogAll() << "process-id: " << chi::mpi.location_id;

  return 0;
}

void chi_lua_test::lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiLuaTest);
}

