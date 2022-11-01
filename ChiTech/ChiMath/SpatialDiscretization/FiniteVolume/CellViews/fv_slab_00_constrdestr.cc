#include "fv_slab.h"

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

/**Slab FV mapping constructor.*/
chi_math::SlabFVValues::
  SlabFVValues(const chi_mesh::Cell& slab_cell,
               const chi_mesh::MeshContinuumConstPtr& grid) :
  CellFVValues(grid,
               slab_cell.centroid,
               std::vector<std::vector<int>>(slab_cell.faces.size(),{-1}))
{
  v0i = slab_cell.vertex_ids[0];
  v1i = slab_cell.vertex_ids[1];
  const auto& v0 = grid->vertices[v0i];
  const auto& v1 = grid->vertices[v1i];

  chi_mesh::Vector3 v01 = v1 - v0;
  volume = v01.Norm();
  face_area.push_back(1.0);
  face_area.push_back(1.0);
}
