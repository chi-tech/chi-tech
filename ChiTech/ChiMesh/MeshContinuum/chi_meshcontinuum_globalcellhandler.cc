#include "chi_meshcontinuum.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

//###################################################################
/**Adds a new cell to grid registry.*/
void chi_mesh::GlobalCellHandler::
  push_back(chi_mesh::Cell *new_cell)
{
  if (new_cell->partition_id == static_cast<uint64_t>(chi_mpi.location_id))
  {
    local_cell_glob_indices.push_back(new_cell->global_id);
    size_t local_cell_index = local_cell_glob_indices.size() - 1;
    new_cell->local_id = local_cell_index;

    native_cells.push_back(new_cell);

    global_cell_id_to_native_id_map.insert(std::make_pair(
      new_cell->global_id, native_cells.size()-1));
  }
  else
  {
    foreign_cells.push_back(new_cell);

    global_cell_id_to_foreign_id_map.insert(std::make_pair(
      new_cell->global_id, foreign_cells.size() - 1));
  }

}

//###################################################################
/**Returns a reference to a cell given its global cell index.*/
chi_mesh::Cell& chi_mesh::GlobalCellHandler::
  operator[](uint64_t cell_global_index)
{
  auto native_location = global_cell_id_to_native_id_map.find(cell_global_index);

  if (native_location != global_cell_id_to_native_id_map.end())
    return *native_cells[native_location->second];
  else
  {
    auto foreign_location = global_cell_id_to_foreign_id_map.find(cell_global_index);
    if (foreign_location != global_cell_id_to_foreign_id_map.end())
      return *foreign_cells[foreign_location->second];
  }

  std::stringstream ostr;
  ostr << "chi_mesh::MeshContinuum::cells. Mapping error."
       << "\n"
       << cell_global_index;

  throw std::invalid_argument(ostr.str());
}

//###################################################################
/**Returns a const reference to a cell given its global cell index.*/
const chi_mesh::Cell& chi_mesh::GlobalCellHandler::
  operator[](uint64_t cell_global_index) const
{
  auto native_location = global_cell_id_to_native_id_map.find(cell_global_index);

  if (native_location != global_cell_id_to_native_id_map.end())
    return *native_cells[native_location->second];
  else
  {
    auto foreign_location = global_cell_id_to_foreign_id_map.find(cell_global_index);
    if (foreign_location != global_cell_id_to_foreign_id_map.end())
      return *foreign_cells[foreign_location->second];
  }

  std::stringstream ostr;
  ostr << "chi_mesh::MeshContinuum::cells. Mapping error."
       << "\n"
       << cell_global_index;

  throw std::invalid_argument(ostr.str());
}

//###################################################################
/**Returns the total number of global cells.*/
size_t chi_mesh::MeshContinuum::GetGlobalNumberOfCells() const
{
  size_t num_local_cells = native_cells.size();
  size_t num_globl_cells = 0;

  MPI_Allreduce(&num_local_cells,
                &num_globl_cells,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD);

  return num_globl_cells;
}

//###################################################################
/**Returns the cell global ids of all ghost cells. These are cells that
 * neighbors to this partition's cells but are on a different
 * partition.*/
std::vector<uint64_t> chi_mesh::GlobalCellHandler::
  GetGhostGlobalIDs()
{
  std::vector<uint64_t> ids;
  ids.reserve(GetNumGhosts());

  for (auto& cell : foreign_cells)
    ids.push_back(cell->global_id);

  return ids;
}

//###################################################################
/**Returns the local storage address of a ghost cell. If the
 * ghost is not truly a ghost then -1 is returned, but is wasteful and
 * therefore the user of this function should implement code
 * to prevent it.*/
uint64_t chi_mesh::GlobalCellHandler::
  GetGhostLocalID(int cell_global_index)
{
  auto foreign_location = global_cell_id_to_foreign_id_map.find(cell_global_index);

  if (foreign_location != global_cell_id_to_foreign_id_map.end())
    return foreign_location->second;

  std::stringstream ostr;
  ostr << "Grid GetGhostLocalID failed to find cell " << cell_global_index;

  throw std::invalid_argument(ostr.str());
}