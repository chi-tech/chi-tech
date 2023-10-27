#include "PieceWiseLinearSlabMapping.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/SpatialDiscretization/CellMappings/PieceWiseLinearBaseMapping.h"

#include "chi_log.h"

namespace chi_math::cell_mapping
{

PieceWiseLinearSlabMapping::PieceWiseLinearSlabMapping(
  const chi_mesh::Cell& slab_cell,
  const chi_mesh::MeshContinuum& ref_grid,
  const chi_math::QuadratureLine& volume_quadrature,
  CoordinateSystemType coordinate_system_type)
  : PieceWiseLinearBaseMapping(ref_grid,
                               slab_cell,
                               2, // num_nodes
                               MakeFaceNodeMapping(slab_cell),
                               coordinate_system_type),
    volume_quadrature_(volume_quadrature)
{
  v0i_ = slab_cell.vertex_ids_[0];
  v1i_ = slab_cell.vertex_ids_[1];
  v0_ = ref_grid_.vertices[v0i_];
  v1_ = ref_grid_.vertices[v1i_];

  chi_mesh::Vector3 v01 = v1_ - v0_;
  h_ = v01.Norm();

  normals_[0] = slab_cell.faces_[0].normal_;
  normals_[1] = slab_cell.faces_[1].normal_;
}

void PieceWiseLinearSlabMapping::ComputeCellVolumeAndAreas(
  const chi_mesh::MeshContinuum& grid,
  const chi_mesh::Cell& cell,
  double& volume,
  std::vector<double>& areas)
{
  std::function<double(const chi_mesh::Vector3&)> swf;
  if (coordinate_system_type_ == CoordinateSystemType::CARTESIAN)
    swf = SpatialDiscretization::CartesianSpatialWeightFunction;
  else if (coordinate_system_type_ == CoordinateSystemType::SPHERICAL)
    swf = SpatialDiscretization::Spherical1DSpatialWeightFunction;
  else
    ChiInvalidArgument("Unsupported coordinate system encountered");

  typedef chi_mesh::Vector3 Vec3;

  volume = 0.0;
  {
    const size_t num_qpoints = volume_quadrature_.qpoints_.size();

    for (size_t qp = 0; qp < num_qpoints; ++qp)
    {
      const Vec3& qpoint = volume_quadrature_.qpoints_[qp];
      const double detJ = h_;

      volume += swf(qpoint) * detJ * volume_quadrature_.weights_[qp];
    }
  } // volume

  areas = {swf(v0_), swf(v1_)};
}

} // namespace chi_math::cell_mapping
