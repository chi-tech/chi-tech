#include "sweep_namespace.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log_exceptions.h"

// ###################################################################
/**Populates the local sub-grid connection information for sweep orderings.*/
void chi_mesh::sweep_management::PopulateCellRelationships(
  const chi_mesh::MeshContinuum& grid,
  const chi_mesh::Vector3& omega,
  std::set<int>& location_dependencies,
  std::set<int>& location_successors,
  std::vector<std::set<std::pair<int, double>>>& cell_successors,
  std::vector<std::vector<FaceOrientation>>& cell_face_orientations)
{
  constexpr double tolerance = 1.0e-16;

  constexpr auto FOPARALLEL = FaceOrientation::PARALLEL;
  constexpr auto FOINCOMING = FaceOrientation::INCOMING;
  constexpr auto FOOUTGOING = FaceOrientation::OUTGOING;

  cell_face_orientations.assign(grid.local_cells.size(), {});
  for (auto& cell : grid.local_cells)
    cell_face_orientations[cell.local_id_].assign(cell.faces_.size(),
                                                  FOPARALLEL);

  for (auto& cell : grid.local_cells)
  {
    size_t f = 0;
    for (auto& face : cell.faces_)
    {
      //======================================= Determine if the face
      //                                        is incident
      FaceOrientation orientation = FOPARALLEL;
      const double mu = omega.Dot(face.normal_);

      bool owns_face = true;
      if (face.has_neighbor_ and cell.global_id_ > face.neighbor_id_ and grid.IsCellLocal(face.neighbor_id_))
        owns_face = false;

      if (owns_face)
      {
        // clang-format off
        if (mu > tolerance) orientation = FOOUTGOING;
        else if (mu < tolerance) orientation = FOINCOMING;

        cell_face_orientations[cell.local_id_][f] = orientation;

        if (face.has_neighbor_ and grid.IsCellLocal(face.neighbor_id_))
        {
          const auto& adj_cell = grid.cells[face.neighbor_id_];
          const auto ass_face = face.GetNeighborAssociatedFace(grid);
          auto& adj_face_ori =
            cell_face_orientations[adj_cell.local_id_][ass_face];

          switch (orientation)
          {
            case FOPARALLEL: adj_face_ori = FOPARALLEL; break;
            case FOINCOMING: adj_face_ori = FOOUTGOING; break;
            case FOOUTGOING: adj_face_ori = FOINCOMING; break;
          }
        }
        // clang-format on
      } // if face owned
      else if (face.has_neighbor_ and not grid.IsCellLocal(face.neighbor_id_))
      {
        const auto& adj_cell = grid.cells[face.neighbor_id_];
        const auto ass_face = face.GetNeighborAssociatedFace(grid);
        const auto& adj_face = adj_cell.faces_[ass_face];

        auto& cur_face_ori =
          cell_face_orientations[cell.local_id_][f];

        const double adj_mu = omega.Dot(adj_face.normal_);
        if (adj_mu > tolerance) orientation = FOOUTGOING;
        else if (adj_mu < tolerance) orientation = FOINCOMING;

        switch (orientation)
        {
          case FOPARALLEL: cur_face_ori = FOPARALLEL; break;
          case FOINCOMING: cur_face_ori = FOOUTGOING; break;
          case FOOUTGOING: cur_face_ori = FOINCOMING; break;
        }
      } // if not face owned locally at all

      ++f;
    } // for face
  }

  //============================================= Make directed connections
  for (auto& cell : grid.local_cells)
  {
    const uint64_t c = cell.local_id_;
    size_t f = 0;
    for (auto& face : cell.faces_)
    {
      const double mu = omega.Dot(face.normal_);
      //======================================= If outgoing determine if
      //                                        it is to a local cell
      if (cell_face_orientations[cell.local_id_][f] == FOOUTGOING)
      {
        //================================ If it is a cell and not bndry
        if (face.has_neighbor_)
        {
          //========================= If it is in the current location
          if (face.IsNeighborLocal(grid))
          {
            double weight = mu * face.ComputeFaceArea(grid);
            cell_successors[c].insert(
              std::make_pair(face.GetNeighborLocalID(grid), weight));
          }
          else
            location_successors.insert(face.GetNeighborPartitionID(grid));
        }
      }
      //======================================= If not outgoing determine
      //                                        what it is dependent on
      else
      {
        //================================if it is a cell and not bndry
        if (face.has_neighbor_ and not face.IsNeighborLocal(grid))
          location_dependencies.insert(face.GetNeighborPartitionID(grid));
      }
      ++f;
    } // for face
  }   // for cell
}
