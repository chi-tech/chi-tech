#ifndef CHITECH_UNIT_TESTS_H
#define CHITECH_UNIT_TESTS_H

/**Declare unit tests here. They should return a bool and the name
 * should start with Test. They should also be passed a verbosity flag.
 * Contrary to the general coding convention the names of the unit test
 * can be Leading_snake_case or Upper_snake_case.*/
namespace chi_unit_tests
{
  bool Test_chi_math(bool verbose);
  bool Test_chi_misc_utils(bool verbose);
  bool Test_chi_data_types(bool verbose);
}

#endif //CHITECH_UNIT_TESTS_H
