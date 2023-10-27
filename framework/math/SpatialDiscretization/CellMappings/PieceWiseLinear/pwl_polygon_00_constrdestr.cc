#include "PieceWiseLinearPolygonMapping.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/SpatialDiscretization/CellMappings/PieceWiseLinearBaseMapping.h"

#include "chi_log.h"

namespace chi_math::cell_mapping
{

// ###################################################################
/** Constructor.*/
PieceWiseLinearPolygonMapping::PieceWiseLinearPolygonMapping(
  const chi_mesh::Cell& poly_cell,
  const chi_mesh::MeshContinuum& ref_grid,
  const chi_math::QuadratureTriangle& volume_quadrature,
  const chi_math::QuadratureLine& surface_quadrature,
  CoordinateSystemType coordinate_system_type)
  : PieceWiseLinearBaseMapping(ref_grid,
                               poly_cell,
                               poly_cell.vertex_ids_.size(), // num_nodes
                               MakeFaceNodeMapping(poly_cell),
                               coordinate_system_type),
    volume_quadrature_(volume_quadrature),
    surface_quadrature_(surface_quadrature)
{
  num_of_subtris_ = static_cast<int>(poly_cell.faces_.size());
  beta_ = 1.0 / num_of_subtris_;

  //=========================================== Get raw vertices
  vc_ = poly_cell.centroid_;

  //=========================================== Calculate legs and determinants
  for (int side = 0; side < num_of_subtris_; side++)
  {
    const chi_mesh::CellFace& face = poly_cell.faces_[side];

    const auto& v0 = ref_grid_.vertices[face.vertex_ids_[0]];
    const auto& v1 = ref_grid_.vertices[face.vertex_ids_[1]];
    chi_mesh::Vertex v2 = vc_;

    chi_mesh::Vector3 sidev01 = v1 - v0;
    chi_mesh::Vector3 sidev02 = v2 - v0;

    double sidedetJ = ((sidev01.x) * (sidev02.y) - (sidev02.x) * (sidev01.y));

    FEside_data2d triangle_data;
    triangle_data.detJ = sidedetJ;
    triangle_data.detJ_surf = sidev01.Norm();

    triangle_data.v_index[0] = face.vertex_ids_[0];
    triangle_data.v_index[1] = face.vertex_ids_[1];

    triangle_data.v0 = v0;

    // Set Jacobian
    triangle_data.J.SetIJ(0, 0, sidev01.x);
    triangle_data.J.SetIJ(1, 0, sidev01.y);
    triangle_data.J.SetIJ(0, 1, sidev02.x);
    triangle_data.J.SetIJ(1, 1, sidev02.y);
    triangle_data.J.SetIJ(2, 2, 0.0);

    // Set Jacobian inverse
    triangle_data.Jinv.SetIJ(0, 0, sidev02.y / sidedetJ);
    triangle_data.Jinv.SetIJ(1, 0, -sidev01.y / sidedetJ);
    triangle_data.Jinv.SetIJ(0, 1, -sidev02.x / sidedetJ);
    triangle_data.Jinv.SetIJ(1, 1, sidev01.x / sidedetJ);
    triangle_data.Jinv.SetIJ(2, 2, 0.0);

    // Set Jacobian-Transpose inverse
    triangle_data.JTinv.SetIJ(0, 0, sidev02.y / sidedetJ);
    triangle_data.JTinv.SetIJ(1, 0, -sidev02.x / sidedetJ);
    triangle_data.JTinv.SetIJ(0, 1, -sidev01.y / sidedetJ);
    triangle_data.JTinv.SetIJ(1, 1, sidev01.x / sidedetJ);
    triangle_data.JTinv.SetIJ(2, 2, 0.0);

    // Set face normal
    triangle_data.normal = face.normal_;

    sides_.push_back(triangle_data);
  }

  //=========================================== Compute node to side mapping
  for (int v = 0; v < poly_cell.vertex_ids_.size(); v++)
  {
    const uint64_t vindex = poly_cell.vertex_ids_[v];
    std::vector<int> side_mapping(num_of_subtris_);
    for (int side = 0; side < num_of_subtris_; side++)
    {
      side_mapping[side] = -1;

      const chi_mesh::CellFace& face = poly_cell.faces_[side];
      if (face.vertex_ids_[0] == vindex) { side_mapping[side] = 0; }
      if (face.vertex_ids_[1] == vindex) { side_mapping[side] = 1; }
    }
    node_to_side_map_.push_back(side_mapping);
  }
}

void PieceWiseLinearPolygonMapping::ComputeCellVolumeAndAreas(
  const chi_mesh::MeshContinuum& grid,
  const chi_mesh::Cell& cell,
  double& volume,
  std::vector<double>& areas)
{
  std::function<double(const chi_mesh::Vector3&)> swf;
  if (coordinate_system_type_ == CoordinateSystemType::CARTESIAN)
    swf = SpatialDiscretization::CartesianSpatialWeightFunction;
  else if (coordinate_system_type_ == CoordinateSystemType::CYLINDRICAL)
    swf = SpatialDiscretization::CylindricalRZSpatialWeightFunction;
  else
    ChiInvalidArgument("Unsupported coordinate system encountered");

  typedef chi_mesh::Vector3 Vec3;

  volume = 0.0;
  for (const auto& side : sides_)
  {
    const size_t num_qpoints = volume_quadrature_.qpoints_.size();

    for (size_t qp = 0; qp < num_qpoints; ++qp)
    {
      const Vec3& qpoint = volume_quadrature_.qpoints_[qp];
      const double detJ = side.detJ;

      const Vec3 qpoint_xyz = side.v0 + side.J * qpoint;

      volume += swf(qpoint_xyz) * detJ * volume_quadrature_.weights_[qp];
    }
  } // for side s

  areas.clear();
  {
    const size_t num_faces = cell.faces_.size();
    for (size_t f = 0; f < num_faces; ++f)
    {
      const size_t num_qpoints = surface_quadrature_.qpoints_.size();

      double area = 0.0;
      for (size_t qp = 0; qp < num_qpoints; ++qp)
      {
        const Vec3& qpoint_face = surface_quadrature_.qpoints_[qp];
        const double detJ = sides_[f].detJ_surf;

        const Vec3 qpoint_xyz = sides_[f].v0 + sides_[f].J * qpoint_face;

        area += swf(qpoint_xyz) * detJ * surface_quadrature_.weights_[qp];
      }
      areas.push_back(area);
    } // for face
  }   // areas
}

} // namespace chi_math::cell_mapping
