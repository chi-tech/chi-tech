#include "chi_lua.h"
#include "chi_runtime.h"

#include "lua_test.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)

#include "unit_tests.h"
#include "chi_log.h"

#include <stdexcept>

//###################################################################
/**This is a lua test function.
\param argument1 Any Argument of any type.
\ingroup LuaGeneralUtilities
 */
int chiLuaTest(lua_State* L)
{
  const int num_args = lua_gettop(L);
  chi::log.Log() << "Hello from chiLuaTest()";
  chi::log.Log() << "num_args = " << num_args;
  return 0;
}

/**Throws an exception that will cause ChiTech to stop
 * execution.*/
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

/**Throws an exception that will not cause ChiTech to stop execution.*/
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

/**Registers the lua function calls for unit tests.*/
void chi_lua_test::lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiLuaTest);
  LUA_FMACRO1(chiThrowException);
  LUA_FMACRO1(chiThrowRecoverableException);

  lua_register(L, "chiUnitTests_Test_chi_math", chi_unit_tests::Test_chi_math);
  lua_register(L, "chiUnitTests_Test_chi_misc_utils",
               chi_unit_tests::Test_chi_misc_utils);
  lua_register(L, "chiUnitTests_Test_chi_data_types",
               chi_unit_tests::Test_chi_data_types);
  lua_register(L, "chiUnitTests_Test_WDD_IJK_Sweep",
               chi_unit_tests::Test_WDD_IJK_Sweep);
  lua_register(L, "chiUnitTests_Test_paramblock",
               chi_unit_tests::Test_paramblock);

  lua_register(L, "chiSimTest01_FV"  , chi_unit_sim_tests::chiSimTest01_FV);
  lua_register(L, "chiSimTest02_FV"  , chi_unit_sim_tests::chiSimTest02_FV);
  lua_register(L, "chiSimTest03_PWLC", chi_unit_sim_tests::chiSimTest03_PWLC);
  lua_register(L, "chiSimTest04_PWLC", chi_unit_sim_tests::chiSimTest04_PWLC);

  lua_register(L, "chiSimTest06_WDD", chi_unit_sim_tests::chiSimTest06_WDD);

  lua_register(L, "chiSimTest91_PWLD", chi_unit_sim_tests::chiSimTest91_PWLD);
  lua_register(L, "chiSimTest92_DSA", chi_unit_sim_tests::chiSimTest92_DSA);
  lua_register(L, "chiSimTest93_RayTracing",
               chi_unit_sim_tests::chiSimTest93_RayTracing);

  lua_register(L, "chiSimTest_IP_MMS_L2error",
               chi_unit_sim_tests::chiSimTest_IP_MMS_L2error);
}

