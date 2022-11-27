#include "pwl_slab.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "pwl_cellbase.h"

chi_math::SlabMappingFE_PWL::
SlabMappingFE_PWL(const chi_mesh::Cell& slab_cell,
                  const chi_mesh::MeshContinuumConstPtr& ref_grid,
                  const chi_math::QuadratureLine &volume_quadrature) :
  chi_math::CellMappingFE_PWL(ref_grid,
                              slab_cell,
                              2, //num_nodes
                              GetVertexLocations(*ref_grid,slab_cell),
                              MakeFaceNodeMapping(slab_cell)),
  volume_quadrature(volume_quadrature)
{
  v0i = slab_cell.vertex_ids[0];
  v1i = slab_cell.vertex_ids[1];
  v0 = m_grid_ptr->vertices[v0i];
  const auto& v1 = m_grid_ptr->vertices[v1i];

  chi_mesh::Vector3 v01 = v1 - v0;
  h = v01.Norm();

//  face_node_mappings.emplace_back(1, 0);
//  face_node_mappings.emplace_back(1, 1);

  normals[0] = slab_cell.faces[0].normal;
  normals[1] = slab_cell.faces[1].normal;
}