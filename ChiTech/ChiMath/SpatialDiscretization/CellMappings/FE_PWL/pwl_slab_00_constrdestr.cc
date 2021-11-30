#include "pwl_slab.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

SlabMappingFE_PWL::
SlabMappingFE_PWL(const chi_mesh::Cell& slab_cell,
                  const chi_mesh::MeshContinuumPtr &ref_grid,
                  const chi_math::QuadratureLine &volume_quadrature) :
  CellMappingFE_PWL(2, ref_grid),
  volume_quadrature(volume_quadrature)
{
  grid = ref_grid;
  v0i = slab_cell.vertex_ids[0];
  v1i = slab_cell.vertex_ids[1];
  v0 = grid->vertices[v0i];
  const auto& v1 = grid->vertices[v1i];

  chi_mesh::Vector3 v01 = v1 - v0;
  h = v01.Norm();

  face_dof_mappings.emplace_back(1,0);
  face_dof_mappings.emplace_back(1,1);

  normals[0] = slab_cell.faces[0].normal;
  normals[1] = slab_cell.faces[1].normal;
}