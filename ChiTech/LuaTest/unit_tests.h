#ifndef CHITECH_UNIT_TESTS_H
#define CHITECH_UNIT_TESTS_H

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

#endif //CHITECH_UNIT_TESTS_H
