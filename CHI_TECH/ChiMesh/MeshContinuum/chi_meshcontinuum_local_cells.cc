#include "chi_meshcontinuum.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

//###################################################################
/**Returns a reference to a local cell, given a local cell index.*/
chi_mesh::Cell& chi_mesh::MeshContinuum::LocalCells::
operator[](int cell_local_index)
{
  if ((cell_local_index < 0) or
      (cell_local_index >= local_cell_ind.size()))
  {
    chi_log.Log(LOG_ALLERROR)
      << "LocalCells attempted to access local cell "
      << cell_local_index
      << " but index out of range [0, "
      << local_cell_ind.size()-1 << "].";
    exit(EXIT_FAILURE);
  }
//  int cell_global_index = local_cell_ind[cell_local_index];
//  return *cell_references[cell_global_index];
  return *native_cells[cell_local_index];
}

