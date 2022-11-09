#include "spatial_discretization.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

const chi_math::CellMapping& chi_math::SpatialDiscretization::
  GetCellMapping(const chi_mesh::Cell& cell) const
{
  constexpr std::string_view fname = "chi_math::SpatialDiscretization::"
                                     "GetCellMapping";
  try
  {
    if (ref_grid->IsCellLocal(cell.global_id))
    {
      return *cell_mappings.at(cell.local_id);
    }
    else
    {
      return *nb_cell_mappings.at(cell.global_id);
    }
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range(
      std::string(fname) + ": Failed to obtain cell mapping.");
  }
}