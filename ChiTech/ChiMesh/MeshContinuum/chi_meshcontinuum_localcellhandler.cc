#include "chi_meshcontinuum.h"

//###################################################################
/**Returns a reference to a local cell, given a local cell index.*/
chi_mesh::Cell& chi_mesh::LocalCellHandler::
  operator[](uint64_t cell_local_index)
{
  if (cell_local_index >= native_cells.size())
  {
    std::stringstream ostr;
    ostr
      << "LocalCells attempted to access local cell "
      << cell_local_index
      << " but index out of range [0, "
      << native_cells.size() - 1 << "].";
    throw std::invalid_argument(ostr.str());
  }

  return *native_cells[cell_local_index];
}

