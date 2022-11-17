#ifndef CHITECH_UNIT_TESTS_H
#define CHITECH_UNIT_TESTS_H

#include "chi_lua.h"

#define ChiUnitTestMessageHome(cond)                              \
{                                                            \
  if (not cond)                                              \
    chi::log.Log0() << __PRETTY_FUNCTION__ << " Failed";   \
  else                                                       \
    chi::log.Log0() << __PRETTY_FUNCTION__ << " Succeeded";\
}
#define ChiUnitTestMessageAll(cond)                              \
{                                                            \
  if (not cond)                                              \
    chi::log.LogAll() << __PRETTY_FUNCTION__ << " Failed";   \
  else                                                       \
    chi::log.LogAll() << __PRETTY_FUNCTION__ << " Succeeded";\
}
/**Declare unit tests here. They should return a bool and the name
 * should start with Test. They should also be passed a verbosity flag.
 * Contrary to the general coding convention the names of the unit test
 * can be Leading_snake_case or Upper_snake_case.*/
namespace chi_unit_tests
{
  void Test_chi_math(bool verbose);
  void Test_chi_misc_utils(bool verbose);
  void Test_chi_data_types(bool verbose);
  void Test_WDD_IJK_Sweep(bool verbose);
}

/**Declare unit tests here that are meant to run from the lua console.*/
namespace chi_unit_sim_tests
{
  int chiSimTest01_FV(lua_State* L);
  int chiSimTest02_FV(lua_State* L);
  int chiSimTest03_PWLC(lua_State* L);
  int chiSimTest04_PWLC(lua_State* L);

  int chiSimTest91_PWLD(lua_State* L);
}
#endif //CHITECH_UNIT_TESTS_H
