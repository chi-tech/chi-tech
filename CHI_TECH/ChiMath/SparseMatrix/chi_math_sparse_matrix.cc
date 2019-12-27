#include "chi_math_sparse_matrix.h"

#include <iomanip>
#include <unistd.h>

#include <chi_log.h>

extern ChiLog chi_log;



//###################################################################
/**Default constructor.*/
chi_math::SparseMatrix::SparseMatrix()
{
  nnz = 0;
}

//###################################################################
/**Copy constructor.*/
chi_math::SparseMatrix::
  SparseMatrix(const chi_math::SparseMatrix& in_matrix)
{
  nnz = in_matrix.nnz;
  row_size = in_matrix.row_size;
  col_size = in_matrix.col_size;

  nnz = 0;
  rowI_colJ.resize(row_size,std::vector<double>());
  inds_rowI.resize(row_size,std::vector<int>());

  j_jr_maps.resize(row_size,std::vector<int>(col_size,-1));

  for (int i=0; i<in_matrix.rowI_colJ.size(); i++)
  {
    rowI_colJ[i] = (in_matrix.rowI_colJ[i]);
    inds_rowI[i] = (in_matrix.inds_rowI[i]);
    j_jr_maps[i] = (in_matrix.j_jr_maps[i]);
  }
}


//###################################################################
/**Constructor with number of rows constructor.*/
chi_math::SparseMatrix::SparseMatrix(int num_rows, int num_cols)
{
  nnz = 0;
  rowI_colJ.resize(num_rows,std::vector<double>());
  inds_rowI.resize(num_rows,std::vector<int>());

  j_jr_maps.resize(num_rows,std::vector<int>(num_cols,-1));

  row_size = num_rows;
  col_size = num_cols;
}

//###################################################################
/**Obtains a reference to an ij value by first mapping the j.*/
double& chi_math::SparseMatrix::GetIJ(int i, int j)
{
  int jr = j_jr_maps[i][j];

  return rowI_colJ[i][jr];
}

//###################################################################
/**Inserts a value into the matrix with no duplicate check.*/
void chi_math::SparseMatrix::SetSize(int num_rows, int num_cols)
{
  nnz = 0;
  rowI_colJ.resize(num_rows,std::vector<double>());
  inds_rowI.resize(num_rows,std::vector<int>());

  j_jr_maps.resize(num_rows,std::vector<int>(num_cols,-1));

  row_size = num_rows;
  col_size = num_cols;
}

//###################################################################
/**Inserts a value into the matrix with no duplicate check.*/
void chi_math::SparseMatrix::UnsafeInsert(int i, int j, double value)
{
  CheckInitialized();

  rowI_colJ[i].push_back(value);
  inds_rowI[i].push_back(j);
  j_jr_maps[i][j] = inds_rowI[i].size()-1;
}

//###################################################################
/**Inserts a value into the matrix with duplicate check.*/
void chi_math::SparseMatrix::Insert(int i, int j, double value)
{
  CheckInitialized();

  if ((i<0) || (i>=row_size) || (j<0) || (j>=col_size))
  {
    chi_log.Log(LOG_ALLERROR)
      << "SparseMatrix::Insert encountered out of bounds,"
      << " i=" << i << " j=" << j
      << " bounds(" << row_size << "," << col_size << ")";
    exit(EXIT_FAILURE);
  }

  bool already_there = false;
  if (j_jr_maps[i][j]>=0)
  {
    already_there = true;
    int jr = j_jr_maps[i][j];
    rowI_colJ[i][jr] = value;
  }

  if (!already_there)
  {
    rowI_colJ[i].push_back(value);
    inds_rowI[i].push_back(j);
    j_jr_maps[i][j] = inds_rowI[i].size()-1;
  }
}

//###################################################################
/**Inserts-Adds a value into the matrix with duplicate check.*/
void chi_math::SparseMatrix::InsertAdd(int i, int j, double value)
{
  CheckInitialized();

  if ((i<0) || (i>=row_size) || (j<0) || (j>=col_size))
  {
    chi_log.Log(LOG_ALLERROR)
      << "SparseMatrix::Insert encountered out of bounds,"
      << " i=" << i << " j=" << j
      << " bounds(" << row_size << "," << col_size << ")";
    exit(EXIT_FAILURE);
  }

  bool already_there = false;
  if (j_jr_maps[i][j]>=0)
  {
    already_there = true;
    int jr = j_jr_maps[i][j];
    rowI_colJ[i][jr] += value;
  }

  if (!already_there)
  {
    rowI_colJ[i].push_back(value);
    inds_rowI[i].push_back(j);
    j_jr_maps[i][j] = inds_rowI[i].size()-1;
  }
}

//###################################################################
/**Sets the diagonal of the matrix using a vector.*/
void chi_math::SparseMatrix::SetDiagonal(const std::vector<double>& diag)
{
  CheckInitialized();

  int num_rows = rowI_colJ.size();
  //============================================= Check size
  if (diag.size() != rowI_colJ.size())
  {
    chi_log.Log(LOG_ALLERROR)
    << "Incompatible matrix-vector size encountered "
    << "in call to SparseMatrix::SetDiagonal.";
    exit(EXIT_FAILURE);
  }

  //============================================= Assign values
  for (int i=0; i<num_rows; i++)
  {
    int num_indices = rowI_colJ[i].size();
    bool already_there = false;
    for (int jr=0; jr<num_indices; jr++)
    {
      if (inds_rowI[i][jr] == i)
      {
        rowI_colJ[i][jr] = diag[i];
        already_there = true;
      }
    }

    if (!already_there)
    {
      rowI_colJ[i].push_back(diag[i]);
      inds_rowI[i].push_back(i);
      j_jr_maps[i][i] = inds_rowI[i].size()-1;
    }
  }
}

//###################################################################
/**Bracket operator for the matrix using helper.*/
chi_math::SparseMatrix::SparseMatrixHelper
  chi_math::SparseMatrix::operator[](int i)
{
  CheckInitialized();

  if ((i <0) || (i>=row_size))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Sparse Matrix row index out of bounds"
      << " in call to SparseMatrix::operator[]";
    exit(EXIT_FAILURE);
  }

  return SparseMatrixHelper(*this,i);
}

//###################################################################
/**Bracket operator for the matrix helper.*/
double& chi_math::SparseMatrix::SparseMatrixHelper::operator[](int j)
{
  if ((j <0) || (j>=parent.col_size))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Sparse Matrix column index out of bounds"
      << " in call to SparseMatrixHelper::operator[]"
      << " j=" << j << " J=" << parent.rowI_colJ[ri].size();
    exit(EXIT_FAILURE);
  }

  if (parent.j_jr_maps[ri][j]<0)
  {
    parent.Insert(ri,j,0.0);

    return parent.GetIJ(ri,j);
  }

  return parent.GetIJ(ri,j);
}

//###################################################################
/**Returns the value in the matrix at the given location.*/
double chi_math::SparseMatrix::ValueIJ(int i, int j)
{
  double retval = 0.0;
  if ((i<0) || (i>=inds_rowI.size()))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Index i out of bounds"
      << " in call to SparseMatrix::ValueIJ"
      << " i=" << i;
    exit(EXIT_FAILURE);
  }

  for (int jr=0; jr<inds_rowI[i].size(); jr++)
  {
    if (inds_rowI[i][jr] == j)
    {
      retval = rowI_colJ[i][jr];
      break;
    }
  }
  return retval;
}

//###################################################################
/**Prints the sparse matrix to string.*/
std::string chi_math::SparseMatrix::PrintS()
{
  std::stringstream out;


  for (int i=0; i<row_size; i++)
  {
    for (int j=0; j<col_size; j++)
    {
      if (j_jr_maps[i][j]<0)
      {
        out
          << std::setprecision(0)
          << std::fixed
          << std::setw(9)
          << 0.0 << " ";
      }
      else
      {
        out
          << std::setprecision(2)
          << std::scientific
          << std::setw(9)
          << rowI_colJ[i][j_jr_maps[i][j]] << " ";
      }
    }
    out << "\n";
  }

  return out.str();
}

//###################################################################
/**Constructor with number of rows constructor.*/
void chi_math::SparseMatrix::CheckInitialized()
{
  if (rowI_colJ.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Illegal call to unitialized SparseMatrix matrix.";
    exit(EXIT_FAILURE);
  }
}
