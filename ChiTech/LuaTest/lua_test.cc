#include <ChiLua/chi_lua.h>

#include "lua_test.h"

#include "unit_tests.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

#define LUA_FMACRO1(x) lua_register(L, #x, x)

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

  if (not chi_unit_tests::Test_chi_math(verbose))
    chi_log.Log(LOG_ALL) << "chi_unit_tests::Test_chi_math Failed";
  if (not chi_unit_tests::Test_chi_misc_utils(verbose))
    chi_log.Log(LOG_ALL) << "chi_unit_tests::Test_chi_misc_utils Failed";
  if (not chi_unit_tests::Test_chi_data_types(verbose))
    chi_log.Log(LOG_ALL) << "chi_unit_tests::Test_chi_data_types Failed";

  return 0;
}

void chi_lua_test::lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiLuaTest);
}

