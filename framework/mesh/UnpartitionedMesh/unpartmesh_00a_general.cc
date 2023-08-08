#include "chi_unpartitioned_mesh.h"

#include "mesh/Cell/cell.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "utils/chi_timer.h"

/**Destructor.*/
chi_mesh::UnpartitionedMesh::~UnpartitionedMesh()
{
  for (auto& cell : raw_cells_)          delete cell;
  for (auto& cell : raw_boundary_cells_) delete cell;
}

/**Compute centroids for all cells.*/
void chi_mesh::UnpartitionedMesh::ComputeCentroidsAndCheckQuality()
{
  const chi_mesh::Vector3 khat(0.0,0.0,1.0);

  Chi::log.Log() << "Computing cell-centroids.";
  for (auto cell : raw_cells_)
  {
    cell->centroid = chi_mesh::Vertex(0.0,0.0,0.0);
    for (auto vid : cell->vertex_ids)
      cell->centroid += vertices_[vid];

    cell->centroid = cell->centroid/static_cast<double>(cell->vertex_ids.size());
  }
  Chi::log.Log() << "Done computing cell-centroids.";

  Chi::log.Log() << "Checking cell-center-to-face orientations";
  size_t num_negative_volume_elements=0;
  for (auto cell : raw_cells_)
  {
    if (cell->type == CellType::POLYGON)
    {
      // Form triangles
      size_t num_verts = cell->vertex_ids.size();
      for (size_t v=0; v<num_verts; ++v)
      {
        size_t vp1 = (v<(num_verts-1))? v+1 : 0;

        const auto& v0 = vertices_[cell->vertex_ids[v]];
        const auto& v1 = vertices_[cell->vertex_ids[vp1]];

        auto E01 = v1 - v0;
        auto n   = E01.Cross(khat).Normalized();

        if (n.Dot(v0 - cell->centroid) < 0.0)
          ++num_negative_volume_elements;
      }//for v
    }
    else if (cell->type == CellType::POLYHEDRON)
    {
      for (auto& face : cell->faces)
      {
        if (face.vertex_ids.size()<2)
          throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                 ": cell-center-to-face check encountered face "
                                 "with less than 2 vertices on a face, making "
                                 "normal computation impossible.");

        // Compute centroid
        chi_mesh::Vector3 face_centroid;
        for (uint64_t vid : face.vertex_ids)
          face_centroid += vertices_[vid];
        face_centroid /= static_cast<double>(face.vertex_ids.size());

        // Form tets for each face edge
        size_t num_face_verts = face.vertex_ids.size();
        for (size_t fv=0; fv<face.vertex_ids.size(); ++fv)
        {
          size_t fvp1 = (fv<(num_face_verts-1))? fv+1 : 0;

          const auto& fv1 = vertices_[face.vertex_ids[fv]];
          const auto& fv2 = vertices_[face.vertex_ids[fvp1]];

          auto E0 = fv1-face_centroid;
          auto E1 = fv2-face_centroid;
          auto n  = E0.Cross(E1).Normalized();

          if (n.Dot(fv1 - cell->centroid) < 0.0)
            ++num_negative_volume_elements;
        }
      }//for face
    }//if polyhedron
  }//for cell in raw_cells

  if (num_negative_volume_elements > 0)
    Chi::log.LogAllWarning()
      << "Cell quality checks detected " << num_negative_volume_elements
      << " negative volume sub-elements (sub-triangle or sub-tetrahedron)."
      << " This issue could result in incorrect quantities"
      << " under some circumstances.";
  Chi::log.Log() << "Done checking cell-center-to-face orientations";
}