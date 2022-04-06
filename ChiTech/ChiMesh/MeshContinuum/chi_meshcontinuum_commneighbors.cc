#include "chi_meshcontinuum.h"

//###################################################################
/**Retrieves a list of all the ghost cells.*/
std::vector<std::unique_ptr<chi_mesh::Cell>> chi_mesh::MeshContinuum::
  GetGhostCells()
{
  const auto& ghost_cell_ids = cells.GetGhostGlobalIDs();

  std::vector<std::unique_ptr<chi_mesh::Cell>> ghost_cells;
  ghost_cells.reserve(ghost_cell_ids.size());

  for (uint64_t global_id : ghost_cell_ids)
  {
    auto cell_copy = std::make_unique<chi_mesh::Cell>(cells[global_id]);
    ghost_cells.push_back(std::move(cell_copy));
  }

  return ghost_cells;
}