#include "pwl_slab.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "pwl_cellbase.h"

chi_math::SlabMappingFE_PWL::
SlabMappingFE_PWL(const chi_mesh::Cell& slab_cell,
                    const chi_mesh::MeshContinuum& ref_grid,
                  const chi_math::QuadratureLine &volume_quadrature) :
    chi_math::CellMappingFE_PWL(ref_grid,
                              slab_cell,
                              2, //num_nodes
                              GetVertexLocations(ref_grid,slab_cell),
                              MakeFaceNodeMapping(slab_cell)),
    volume_quadrature_(volume_quadrature)
{
  v0i_ = slab_cell.vertex_ids_[0];
  v1i_ = slab_cell.vertex_ids_[1];
  v0_ = ref_grid_.vertices[v0i_];
  const auto& v1 = ref_grid_.vertices[v1i_];

  chi_mesh::Vector3 v01 = v1 - v0_;
  h_ = v01.Norm();

//  face_node_mappings.emplace_back(1, 0);
//  face_node_mappings.emplace_back(1, 1);

  normals_[0] = slab_cell.faces_[0].normal_;
  normals_[1] = slab_cell.faces_[1].normal_;
}