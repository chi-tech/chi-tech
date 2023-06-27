#include "spatial_discretization.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

const chi_math::CellMapping& chi_math::SpatialDiscretization::GetCellMapping(
  const chi_mesh::Cell& cell) const
{
  constexpr std::string_view fname = "chi_math::SpatialDiscretization::"
                                     "GetCellMapping";
  try
  {
    if (Grid().IsCellLocal(cell.global_id_))
    {
      return *cell_mappings_.at(cell.local_id_);
    }
    else { return *nb_cell_mappings_.at(cell.global_id_); }
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range(std::string(fname) +
                            ": Failed to obtain cell mapping.");
  }
}

chi_math::SpatialDiscretizationType
chi_math::SpatialDiscretization::Type() const
{
  return type_;
}

const chi_mesh::MeshContinuum& chi_math::SpatialDiscretization::Grid() const
{
  return ref_grid_;
}

chi_math::CoordinateSystemType
chi_math::SpatialDiscretization::GetCoordinateSystemType() const
{
  return coord_sys_type_;
}

size_t chi_math::SpatialDiscretization::GetNumLocalAndGhostDOFs(
  const UnknownManager& unknown_manager) const
{
  return GetNumLocalDOFs(unknown_manager) + GetNumGhostDOFs(unknown_manager);
}