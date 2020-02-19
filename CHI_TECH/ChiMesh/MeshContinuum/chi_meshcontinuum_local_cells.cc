#include "chi_meshcontinuum.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

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
  int cell_global_index = local_cell_ind[cell_local_index];
  return *cell_references[cell_global_index];
}

//###################################################################
/**Adds a new cell to grid registry.*/
void chi_mesh::MeshContinuum::GlobalCellHandler::
  push_back(chi_mesh::Cell *new_cell)
{
  local_cells.cell_references.push_back(new_cell);

  if (new_cell->partition_id == chi_mpi.location_id)
  {
    local_cells.local_cell_ind.push_back(new_cell->cell_global_id);
    int local_cell_index = local_cells.local_cell_ind.size()-1;
    new_cell->cell_local_id = local_cell_index;


    local_cells.native_cells.push_back(new_cell);

    std::pair<int,int> new_info = std::make_pair(
      new_cell->cell_global_id,
      local_cells.native_cells.size()-1);

    global_cell_native_index_set.insert(new_info);
  }
  else
  {
    local_cells.foreign_cells.push_back(new_cell);

    std::pair<int,int> new_info = std::make_pair(
      new_cell->cell_global_id,
      local_cells.foreign_cells.size() - 1);

    global_cell_foreign_index_set.insert(new_info);
  }

}

//###################################################################
/**Returns a pointer to a cell given its global cell index.*/
chi_mesh::Cell* &chi_mesh::MeshContinuum::GlobalCellHandler::
  operator[](int cell_global_index)
{
  return local_cells.cell_references[cell_global_index];
}

//###################################################################
/**Returns the total number of cells.*/
size_t chi_mesh::MeshContinuum::GlobalCellHandler::size()
{
  return local_cells.cell_references.size();
}