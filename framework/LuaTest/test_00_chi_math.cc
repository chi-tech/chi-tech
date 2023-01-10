#include "unit_tests.h"

#include "ChiMath/dynamic_vector.h"
#include "ChiMath/dynamic_matrix.h"
#include "ChiMath/SparseMatrix/chi_math_sparse_matrix.h"

#include "chi_runtime.h"
#include "chi_log.h"

void chi_unit_tests::Test_chi_math(bool verbose)
{
  bool passed = true;
  std::stringstream output;

  //======================================================= Dynamic Vector
  {
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
    } else
      output << std::string("chi_math::DynamicVector<double>.PrintStr() ... Passed\n");
  }
  //======================================================= Dynamic Matrix
  {
    output << "Testing chi_math::DynamicMatrix\n";
    chi_math::DynamicMatrix<double> mat(5, 7, 1.0);

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
    } else
      output << std::string("chi_math::DynamicMatrix<double>.PrintStr() ... Passed\n");
  }

  //======================================================= SparseMatrix
  {
    chi_math::SparseMatrix matrix(4,4);
    auto& mat = matrix;
    mat.Insert(0, 0, 1.0); mat.Insert(0, 1, 1.1); mat.Insert(0, 2, 1.2); mat.Insert(0, 3, 1.3);
    mat.Insert(1, 0, 1.9); mat.Insert(1, 1, 2.0); mat.Insert(1, 2, 2.1);
    mat.Insert(2, 1, 2.9); mat.Insert(2, 2, 3.0); mat.Insert(2, 3, 3.1);
    mat.Insert(3, 2, 3.9); mat.Insert(3, 3, 4.0);

    {
      auto &m = mat;
      output << "----- SparseMatrix::PrintS() -----"
             << "\n" << m.PrintStr() << "\n";

      output << "----- for (const auto& entry : m.Row(2)) -----";
      for (const auto &entry: m.Row(2))
        output << entry.row_index << " "
               << entry.column_index << " "
               << entry.value;

      output << "----- after value*2 -----";
      for (const auto&[row_index, column_index, value]: m.Row(2))
        value *= 2;

      for (const auto &entry: m.Row(2))
        output << entry.row_index << " "
               << entry.column_index << " "
               << entry.value;
    }

    output << "----- for (auto entry : matrix) -----";
    for (const auto& entry : matrix)
      output << entry.row_index << " "
             << entry.column_index << " "
             << entry.value;

    matrix.Compress();
  }

  if (verbose)
    chi::log.Log() << output.str();

  ChiUnitTestMessageHome(passed)
}