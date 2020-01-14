#ifndef _chi_math_sparse_matrix_h
#define _chi_math_sparse_matrix_h

#include "../chi_math.h"

//###################################################################
/**Sparse matrix utility. This is a basic CSR type sparse matrix
 * which allows efficient matrix storage and multiplication. It is
 * not intended for solving linear systems (use PETSc for that instead).
 * It was originally developed for the transfer matrices of transport
 * cross-sections.*/
class chi_math::SparseMatrix
{
private:
  size_t row_size;   ///< Maximum number of rows for this matrix
  size_t col_size;   ///< Maximum number of columns for this matrix

public:
  /**rowI_indices[i] is a vector indices j for the
   * non-zero columns.*/
  std::vector<std::vector<size_t>> rowI_indices;
  /**rowI_values[i] corresponds to column indices and
   * contains the non-zero value.*/
  std::vector<std::vector<double>> rowI_values;

public:
  SparseMatrix(size_t num_rows, size_t num_cols);
  SparseMatrix(const SparseMatrix& in_matrix);

  size_t NumRows() const {return row_size;}
  size_t NumCols() const {return col_size;}

  void   Insert(size_t i, size_t j, double value);
  void   InsertAdd(size_t i, size_t j, double value);
  double ValueIJ(size_t i, size_t j);
  void   SetDiagonal(const std::vector<double>& diag);

  void Compress();

  std::string PrintS();

private:
  void CheckInitialized();
};


#endif