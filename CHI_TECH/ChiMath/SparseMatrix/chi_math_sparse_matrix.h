#ifndef _chi_math_sparse_matrix_h
#define _chi_math_sparse_matrix_h

#include "../chi_math.h"




//###################################################################
/**Sparse matrix utility.*/
class chi_math::SparseMatrix
{
public:
  std::vector<std::vector<double>> rowI_colJ;
  std::vector<std::vector<int>>    inds_rowI;
  std::vector<std::vector<int>>    j_jr_maps;
  double ValueIJ(int i, int j);
private:
  int nnz;
  int row_size;
  int col_size;

  struct SparseMatrixHelper
  {
    chi_math::SparseMatrix& parent;
    int ri;

    SparseMatrixHelper(SparseMatrix& in_parent,int in_ri):parent(in_parent)
    {
      parent = in_parent;
      ri=in_ri;
    }

    double& operator[](int j);

  };

  double& GetIJ(int i,int j);

public:
  SparseMatrix();
  SparseMatrix(const SparseMatrix& in_matrix);
  SparseMatrix(int num_rows, int num_cols);

  void SetSize(int num_rows, int num_cols);
  void UnsafeInsert(int i, int j, double value);
  void Insert(int i, int j, double value);
  void InsertAdd(int i, int j, double value);
//  void SafeInsert(int i,int j, double value);
  void SetDiagonal(const std::vector<double>& diag);

  SparseMatrixHelper operator[](int i);

  std::string PrintS();

private:
  void CheckInitialized();
};


#endif