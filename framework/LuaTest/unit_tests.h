#ifndef CHITECH_UNIT_TESTS_H
#define CHITECH_UNIT_TESTS_H

#include "chi_lua.h"

/**Declare unit tests here.
 * Contrary to the general coding convention the names of the unit test
 * can be Leading_snake_case or Upper_snake_case.*/
namespace chi_unit_tests
{
  int Test_chi_math(lua_State* L);
  int Test_chi_misc_utils(lua_State* L);
  int Test_chi_data_types(lua_State* L);
  int Test_WDD_IJK_Sweep(lua_State* L);
  int Test_paramblock(lua_State* L);
}

/**Declare unit tests here that are meant to run from the lua console.*/
namespace chi_unit_sim_tests
{
  int chiSimTest01_FV(lua_State* L);
  int chiSimTest02_FV(lua_State* L);
  int chiSimTest03_PWLC(lua_State* L);
  int chiSimTest04_PWLC(lua_State* L);

  int chiSimTest06_WDD(lua_State* L);

  int chiSimTest91_PWLD(lua_State* L);
  int chiSimTest92_DSA(lua_State* L);
  int chiSimTest92b_DSA_PWLC(lua_State* L);
  int chiSimTest93_RayTracing(lua_State* L);

  int chiSimTest_IP_MMS_L2error(lua_State* L);
}
#endif //CHITECH_UNIT_TESTS_H
