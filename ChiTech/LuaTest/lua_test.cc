#include "chi_lua.h"
#include "chi_runtime.h"

#include "lua_test.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)

#include "unit_tests.h"

#include <stdexcept>

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
  chi_unit_tests::Test_WDD_IJK_Sweep(verbose);

  return 0;
}


int chiThrowException(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  std::string message = "Unknown exception thrown by " + fname;
  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    message = lua_tostring(L,1);
  }

  throw std::logic_error(message);
}
int chiThrowRecoverableException(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  std::string message = "Unknown exception thrown by " + fname;
  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    message = lua_tostring(L,1);
  }

  throw chi::RecoverableException(message);
}

void chi_lua_test::lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiLuaTest);
  LUA_FMACRO1(chiThrowException);
  LUA_FMACRO1(chiThrowRecoverableException);

  lua_register(L, "chiSimTest01_FV"  , chi_unit_sim_tests::chiSimTest01_FV);
  lua_register(L, "chiSimTest02_FV"  , chi_unit_sim_tests::chiSimTest02_FV);
  lua_register(L, "chiSimTest03_PWLC", chi_unit_sim_tests::chiSimTest03_PWLC);
  lua_register(L, "chiSimTest04_PWLC", chi_unit_sim_tests::chiSimTest04_PWLC);

  lua_register(L, "chiSimTest06_WDD", chi_unit_sim_tests::chiSimTest06_WDD);

  lua_register(L, "chiSimTest91_PWLD", chi_unit_sim_tests::chiSimTest91_PWLD);
}

