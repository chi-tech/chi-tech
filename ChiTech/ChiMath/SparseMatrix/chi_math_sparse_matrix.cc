#include "chi_math_sparse_matrix.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <iomanip>
#include <algorithm>

//###################################################################
/**Constructor with number of rows and columns constructor.*/
chi_math::SparseMatrix::SparseMatrix(size_t num_rows, size_t num_cols) :
  row_size(num_rows),
  col_size(num_cols)
{
  rowI_values.resize(num_rows, std::vector<double>());
  rowI_indices.resize(num_rows, std::vector<size_t>());
}

//###################################################################
/**Copy constructor.*/
chi_math::SparseMatrix::
  SparseMatrix(const chi_math::SparseMatrix& in_matrix)
{
  row_size = in_matrix.NumRows();
  col_size = in_matrix.NumCols();

  rowI_values.resize(row_size, std::vector<double>());
  rowI_indices.resize(row_size, std::vector<size_t>());

  for (size_t i=0; i<in_matrix.rowI_values.size(); i++)
  {
    rowI_values[i] = (in_matrix.rowI_values[i]);
    rowI_indices[i] = (in_matrix.rowI_indices[i]);
  }

}

//###################################################################
/**Inserts a value into the matrix.*/
void chi_math::SparseMatrix::Insert(size_t i, size_t j, double value)
{
  CheckInitialized();

  if ((i<0) || (i>=row_size) || (j<0) || (j>=col_size))
  {
    chi::log.LogAllError()
      << "SparseMatrix::Insert encountered out of bounds,"
      << " i=" << i << " j=" << j
      << " bounds(" << row_size << "," << col_size << ")";
    exit(EXIT_FAILURE);
  }

  auto relative_location = std::find(rowI_indices[i].begin(),
                                     rowI_indices[i].end(), j);
  bool already_there = (relative_location != rowI_indices[i].end());

  if (already_there)
  {
    size_t jr = relative_location - rowI_indices[i].begin();
    rowI_values[i][jr] = value;
  }
  else
  {
    rowI_values[i].push_back(value);
    rowI_indices[i].push_back(j);
  }
}

//###################################################################
/**Inserts-Adds a value into the matrix with duplicate check.*/
void chi_math::SparseMatrix::InsertAdd(size_t i, size_t j, double value)
{
  CheckInitialized();

  if ((i<0) || (i>=row_size) || (j<0) || (j>=col_size))
  {
    chi::log.LogAllError()
      << "SparseMatrix::Insert encountered out of bounds,"
      << " i=" << i << " j=" << j
      << " bounds(" << row_size << "," << col_size << ")";
    exit(EXIT_FAILURE);
  }

  auto relative_location = std::find(rowI_indices[i].begin(),
                                     rowI_indices[i].end(), j);
  bool already_there = (relative_location != rowI_indices[i].end());

  if (already_there)
  {
    size_t jr = relative_location - rowI_indices[i].begin();
    rowI_values[i][jr] += value;
  }
  else
  {
    rowI_values[i].push_back(value);
    rowI_indices[i].push_back(j);
  }
}

//###################################################################
/**Sets the diagonal of the matrix using a vector.*/
void chi_math::SparseMatrix::SetDiagonal(const std::vector<double>& diag)
{
  CheckInitialized();

  size_t num_rows = rowI_values.size();
  //============================================= Check size
  if (diag.size() != rowI_values.size())
  {
    chi::log.LogAllError()
    << "Incompatible matrix-vector size encountered "
    << "in call to SparseMatrix::SetDiagonal.";
    exit(EXIT_FAILURE);
  }

  //============================================= Assign values
  for (size_t i=0; i<num_rows; i++)
  {
    auto relative_location = std::find(rowI_indices[i].begin(),
                                       rowI_indices[i].end(), i);
    bool already_there = (relative_location != rowI_indices[i].end());

    if (already_there)
    {
      size_t jr = relative_location - rowI_indices[i].begin();
      rowI_values[i][jr] = diag[i];
    }
    else
    {
      rowI_values[i].push_back(diag[i]);
      rowI_indices[i].push_back(i);
    }
  }//for i
}

//###################################################################
/**Returns the value in the matrix at the given location. This
 * is a rather inefficient routine. Use the columns and values
 * rather than directly this function.*/
double chi_math::SparseMatrix::ValueIJ(size_t i, size_t j) const
{
  double retval = 0.0;
  if ((i<0) || (i >= rowI_indices.size()))
  {
    chi::log.LogAllError()
      << "Index i out of bounds"
      << " in call to SparseMatrix::ValueIJ"
      << " i=" << i;
    exit(EXIT_FAILURE);
  }

  if (not rowI_indices[i].empty())
  {
    auto relative_location = std::find(rowI_indices[i].begin(),
                                       rowI_indices[i].end(), j);
    bool non_zero = (relative_location != rowI_indices[i].end());
    if (non_zero)
    {
      size_t jr = relative_location - rowI_indices[i].begin();
      retval = rowI_values[i][jr];
    }
  }
  return retval;
}

//###################################################################
/**Sorts the column indices of each row for faster lookup.*/
void chi_math::SparseMatrix::Compress()
{
  for (size_t i=0; i < rowI_indices.size(); ++i)
  {
    auto& indices = rowI_indices[i];
    auto& values  = rowI_values[i];

    //====================================== Copy row indexes and values into
    //                                       vector of pairs
    std::vector<std::pair<size_t,double>> target;
    target.reserve(indices.size());

    auto index = indices.begin();
    auto value = values.begin();
    for (;index!=indices.end(); ++index, ++value)
    {
      target.emplace_back(*index,*value);
    }

    //====================================== Define compare operator
    struct
    {
      bool operator()(std::pair<size_t,double> a, std::pair<size_t,double> b)
      {
        return a.first < b.first;
      }
    }compare_index;

    //Sort
    std::stable_sort(target.begin(), target.end(), compare_index);

    //====================================== Copy back
    indices.clear();
    values.clear();
    for (auto& iv_pair : target)
    {
      indices.push_back(iv_pair.first);
      values.push_back(iv_pair.second);
    }
  }

}

//###################################################################
/**Prints the sparse matrix to string.*/
std::string chi_math::SparseMatrix::PrintS()
{
  std::stringstream out;


  for (size_t i=0; i<row_size; i++)
  {
    for (size_t j=0; j<col_size; j++)
    {
      auto relative_location = std::find(rowI_indices[i].begin(),
                                         rowI_indices[i].end(), j);
      bool non_zero = (relative_location != rowI_indices[i].end());

      if (non_zero)
      {
        size_t jr = relative_location - rowI_indices[i].begin();
        out
          << std::setprecision(2)
          << std::scientific
          << std::setw(9)
          << rowI_values[i][jr] << " ";
      }
      else
      {
        out
          << std::setprecision(0)
          << std::fixed
          << std::setw(9)
          << 0.0 << " ";
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
  if (rowI_values.empty())
  {
    chi::log.LogAllError()
      << "Illegal call to unitialized SparseMatrix matrix.";
    exit(EXIT_FAILURE);
  }
}
