#include "unit_tests.h"

#include "ChiMath/dynamic_vector.h"
#include "ChiMath/dynamic_matrix.h"

#include "chi_log.h"
extern ChiLog& chi_log;

bool chi_unit_tests::Test_chi_math(bool verbose)
{
  bool passed = true;
  std::stringstream output;

  //======================================================= Dynamic Vector
  output << "Testing chi_math::DynamicVector\n";

  chi_math::DynamicVector<double> vec(5, 1.0);

  const auto vec_str = vec.PrintStr();

  if (vec_str != "[1 1 1 1 1]")
  {
    passed = false;
    output << std::string("chi_math::DynamicVector<double>.PrintStr() ... Failed\n"
                          " Expected:\n"
                          "[1 1 1 1 1]\n"
                          "Got:\n") + vec_str;
  }
  else
    output << std::string("chi_math::DynamicVector<double>.PrintStr() ... Passed\n");

  //======================================================= Dynamic Matrix
  output << "Testing chi_math::DynamicMatrix\n";
  chi_math::DynamicMatrix<double> mat(5,7,1.0);

  const auto mat_str = mat.PrintStr();

  if (mat_str != "1 1 1 1 1 1 1\n"
                 "1 1 1 1 1 1 1\n"
                 "1 1 1 1 1 1 1\n"
                 "1 1 1 1 1 1 1\n"
                 "1 1 1 1 1 1 1")
  {
    passed = false;
    output << std::string("chi_math::DynamicMatrix<double>.PrintStr() ... Failed\n"
                          " Expected:\n"
                          "1 1 1 1 1 1 1\n"
                          "1 1 1 1 1 1 1\n"
                          "1 1 1 1 1 1 1\n"
                          "1 1 1 1 1 1 1\n"
                          "1 1 1 1 1 1 1\n"
                          "Got:\n") + mat_str;
  }
  else
    output << std::string("chi_math::DynamicMatrix<double>.PrintStr() ... Passed\n");

  if (verbose)
    chi_log.Log() << output.str();

  return passed;
}