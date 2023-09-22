#include "spatial_discretization.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

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

std::pair<std::set<uint32_t>, std::set<uint32_t>>
chi_math::SpatialDiscretization::MakeCellInternalAndBndryNodeIDs(
  const chi_mesh::Cell& cell) const
{
  const auto& cell_mapping = GetCellMapping(cell);
  const size_t num_faces = cell.faces_.size();
  const size_t num_nodes = cell_mapping.NumNodes();

  //====================================== Determine which nodes are on the
  //                                       boundary
  std::set<uint32_t> boundary_nodes;
  for (size_t f = 0; f < num_faces; ++f)
  {
    if (not cell.faces_[f].has_neighbor_)
    {
      const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
      for (size_t fi = 0; fi < num_face_nodes; ++fi)
        boundary_nodes.insert(cell_mapping.MapFaceNode(f, fi));
    }
  } // for f

  //====================================== Determine non-boundary nodes
  std::set<uint32_t> internal_nodes;
  for (size_t i = 0; i < num_nodes; ++i)
    if (boundary_nodes.find(i) == boundary_nodes.end())
      internal_nodes.insert(i);

  return {internal_nodes, boundary_nodes};
}