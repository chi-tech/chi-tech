#include <ChiLua/chi_lua.h>

#include "lua_test.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)

#include "unit_tests.h"

//###################################################################
/**This is a lua test function.
\param argument1 Any Argument of any type.
\ingroup LuaGeneralUtilities
 */
int chiLuaTest(lua_State* L)
{
  const int num_args = lua_gettop(L);
  bool verbose = false;
  if (num_args >= 1)
    verbose = lua_toboolean(L,1);

  chi_unit_tests::Test_chi_math(verbose);
  chi_unit_tests::Test_chi_misc_utils(verbose);
  chi_unit_tests::Test_chi_data_types(verbose);

  return 0;
}

void chi_lua_test::lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiLuaTest);
}

