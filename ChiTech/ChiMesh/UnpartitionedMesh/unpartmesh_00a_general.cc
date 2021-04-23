#include "chi_unpartitioned_mesh.h"

#include "ChiMesh/Cell/cell_polyhedron.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

/**Compute centroids for all cells.*/
void chi_mesh::UnpartitionedMesh::ComputeCentroidsAndCheckQuality()
{
  chi_log.Log() << "Computing cell-centroids.";
  for (auto cell : raw_cells)
  {
    cell->centroid = chi_mesh::Vertex(0.0,0.0,0.0);
    for (auto vid : cell->vertex_ids)
      cell->centroid += vertices[vid];

    cell->centroid = cell->centroid/double(cell->vertex_ids.size());
  }
  chi_log.Log() << "Done computing cell-centroids.";

  chi_log.Log() << "Checking cell-center-to-face orientations";
  size_t check0_corrections=0;
  for (auto cell : raw_cells)
    if (cell->type == CellType::POLYHEDRON)
      for (auto& face : cell->faces)
      {
        chi_mesh::Vector3 face_centroid;
        for (uint64_t vid : face.vertex_ids)
          face_centroid += vertices[vid];
        face_centroid /= double(face.vertex_ids.size());

        if (face.vertex_ids.size()<2)
          throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                 ": cell-center-to-face check encountered face with less than"
                                 " 2 vertices on a face, making normal computation impossible.");

        const auto& fv1 = vertices[face.vertex_ids[0]];
        const auto& fv2 = vertices[face.vertex_ids[1]];

        auto E0 = fv1-face_centroid;
        auto E1 = fv2-face_centroid;
        auto n  = E0.Cross(E1); n.Normalize();

        auto CC = face_centroid - cell->centroid;

        if (n.Dot(CC) < 0.0)
        {
          std::vector<uint64_t> reversed_verts(face.vertex_ids.rbegin(),
                                               face.vertex_ids.rend());
          face.vertex_ids = reversed_verts;
          ++check0_corrections;
        }
      }//for face

  if (check0_corrections > 0)
    chi_log.Log(LOG_ALLWARNING)
      << "Cell-center-to-face orientation detected " << check0_corrections
      << " faces that violated quality requirements. An attempt to fix it"
      << " was made.";
  chi_log.Log() << "Done checking cell-center-to-face orientations";
}