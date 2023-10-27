#include "LagrangeBaseMapping.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_math::cell_mapping
{

LagrangeBaseMapping::LagrangeBaseMapping(
  const chi_mesh::MeshContinuum& grid,
  const chi_mesh::Cell& cell,
  size_t num_nodes,
  std::vector<std::vector<int>> face_node_mappings,
  const Quadrature& volume_quadrature,
  const Quadrature& surface_quadrature,
  CoordinateSystemType coordinate_system_type)
  : CellMapping(grid,
                cell,
                num_nodes,
                GetVertexLocations(grid, cell),
                std::move(face_node_mappings),
                coordinate_system_type),
    volume_quadrature_(volume_quadrature),
    surface_quadrature_(surface_quadrature)
{
}

/** This section just determines a mapping of face dofs
to cell dofs. This is pretty simple since we can
just loop over each face dof then subsequently
loop over cell dofs, if the face dof node index equals
the cell dof node index then the mapping is assigned.

This mapping is not used by any of the methods in
    this class but is used by methods requiring the
      surface integrals of the shape functions.*/
std::vector<std::vector<int>>
LagrangeBaseMapping::MakeFaceNodeMapping(const chi_mesh::Cell& cell)
{
  const size_t num_faces = cell.faces_.size();
  std::vector<std::vector<int>> mappings;
  mappings.reserve(num_faces);
  for (auto& face : cell.faces_)
  {
    std::vector<int> face_dof_mapping;
    face_dof_mapping.reserve(face.vertex_ids_.size());
    for (uint64_t fvid : face.vertex_ids_)
    {
      int mapping = -1;
      for (size_t ci = 0; ci < cell.vertex_ids_.size(); ci++)
      {
        if (fvid == cell.vertex_ids_[ci])
        {
          mapping = static_cast<int>(ci);
          break;
        }
      } // for cell i
      if (mapping < 0)
      {
        Chi::log.LogAllError() << "Unknown face mapping encountered. "
                                  "pwl_polyhedron.h";
        Chi::Exit(EXIT_FAILURE);
      }
      face_dof_mapping.push_back(mapping);
    } // for face i

    mappings.push_back(face_dof_mapping);
  }
  return mappings;
}

std::vector<chi_mesh::Vector3>
LagrangeBaseMapping::GetVertexLocations(const chi_mesh::MeshContinuum& grid,
                                        const chi_mesh::Cell& cell)
{
  std::vector<chi_mesh::Vector3> verts;
  verts.reserve(cell.vertex_ids_.size());

  for (const auto vid : cell.vertex_ids_)
    verts.push_back(grid.vertices[vid]);

  return verts;
}

chi_mesh::Vector3
LagrangeBaseMapping::MapWorldXYZToNaturalXYZ(const Vec3& world_xyz) const
{
  WorldXYZToNaturalMappingHelper helper(*this, world_xyz);

  auto nat_xyz = chi_math::NewtonIteration(helper, {0.0, 0.0, 0.0}, 20, 1.0e-8);

  return {nat_xyz[0], nat_xyz[1], nat_xyz[2]};
}

double LagrangeBaseMapping::ShapeValue(int i,
                                       const chi_mesh::Vector3& xyz) const
{
  const auto natural_xyz = MapWorldXYZToNaturalXYZ(xyz);
  return RefShape(i, natural_xyz);
}

void LagrangeBaseMapping::ShapeValues(const chi_mesh::Vector3& xyz,
                                      std::vector<double>& shape_values) const
{
  const auto natural_xyz = MapWorldXYZToNaturalXYZ(xyz);

  shape_values.clear();
  for (int i = 0; i < num_nodes_; ++i)
    shape_values.push_back(RefShape(i, natural_xyz));
}

chi_mesh::Vector3
LagrangeBaseMapping::GradShapeValue(int i, const chi_mesh::Vector3& xyz) const
{
  const auto natural_xyz = MapWorldXYZToNaturalXYZ(xyz);
  return RefGradShape(i, natural_xyz);
}

void LagrangeBaseMapping::GradShapeValues(
  const chi_mesh::Vector3& xyz,
  std::vector<chi_mesh::Vector3>& gradshape_values) const
{
  const auto natural_xyz = MapWorldXYZToNaturalXYZ(xyz);

  gradshape_values.clear();
  for (int i = 0; i < num_nodes_; ++i)
    gradshape_values.push_back(RefGradShape(i, natural_xyz));
}

finite_element::VolumetricQuadraturePointData
LagrangeBaseMapping::MakeVolumetricQuadraturePointData() const
{
  typedef std::vector<Vec3> VecVec3;
  const size_t num_qpoints = volume_quadrature_.qpoints_.size();

  //=================================== Declare necessary vars
  std::vector<unsigned int> quadrature_point_indices;
  VecVec3 qpoints_xyz;
  std::vector<VecDbl> shape_value;
  std::vector<VecVec3> shape_grad;
  VecDbl JxW;

  quadrature_point_indices.reserve(num_qpoints);
  for (unsigned int qp = 0; qp < num_qpoints; ++qp)
    quadrature_point_indices.push_back(qp);

  JxW.assign(num_qpoints, 0.0);
  shape_value.assign(num_nodes_, VecDbl(num_qpoints, 0.0));
  shape_grad.assign(num_nodes_, VecVec3(num_qpoints));
  qpoints_xyz.assign(num_qpoints, Vec3{0.0, 0.0, 0.0});

  for (uint32_t qp : quadrature_point_indices)
  {
    const Vec3& qpoint = volume_quadrature_.qpoints_[qp];
    const auto J = RefJacobian(qpoint);
    const auto JT = chi_math::Transpose(J);
    const auto JTinv = chi_math::Inverse(JT);
    const double detJ = Determinant(J);

    JxW[qp] = detJ * volume_quadrature_.weights_[qp];

    for (size_t i = 0; i < num_nodes_; ++i)
    {
      const double ref_shape_i = RefShape(i, qpoint);
      shape_value[i][qp] = ref_shape_i;
      const Vec3 ref_shape_grad = RefGradShape(i, qpoint);
      const VecDbl b = {ref_shape_grad.x, ref_shape_grad.y, ref_shape_grad.z};
      const VecDbl JTInv_b = MatMul(JTinv, b);

      shape_grad[i][qp] = {JTInv_b[0], JTInv_b[1], JTInv_b[2]};

      const auto& node_xyz = node_locations_[i];

      qpoints_xyz[qp] += ref_shape_i * node_xyz;
    } // for i
  }   // for qp

  return finite_element::VolumetricQuadraturePointData(quadrature_point_indices,
                                                       qpoints_xyz,
                                                       shape_value,
                                                       shape_grad,
                                                       JxW,
                                                       face_node_mappings_,
                                                       num_nodes_);
}

const Quadrature&
LagrangeBaseMapping::GetSurfaceQuadrature(size_t face_index) const
{
  return surface_quadrature_;
}

finite_element::SurfaceQuadraturePointData
LagrangeBaseMapping::MakeSurfaceQuadraturePointData(size_t face_index) const
{
  const auto& surface_quadrature = GetSurfaceQuadrature(face_index);

  typedef std::vector<Vec3> VecVec3;
  const size_t num_qpoints = surface_quadrature.qpoints_.size();

  //=================================== Declare necessary vars
  std::vector<unsigned int> quadrature_point_indices;
  VecVec3 qpoints_xyz;
  std::vector<VecDbl> shape_value;
  std::vector<VecVec3> shape_grad;
  VecDbl JxW;
  VecVec3 normals;

  quadrature_point_indices.reserve(num_qpoints);
  for (unsigned int qp = 0; qp < num_qpoints; ++qp)
    quadrature_point_indices.push_back(qp);

  JxW.assign(num_qpoints, 0.0);
  normals.assign(num_qpoints, Vec3{0.0, 0.0, 0.0});
  shape_value.assign(num_nodes_, VecDbl(num_qpoints, 0.0));
  shape_grad.assign(num_nodes_, VecVec3(num_qpoints));
  qpoints_xyz.assign(num_qpoints, Vec3{0.0, 0.0, 0.0});

  const size_t f = face_index;
  for (uint32_t qp : quadrature_point_indices)
  {
    const Vec3& qpoint_face = surface_quadrature.qpoints_[qp];
    const auto qpoint = FaceToElementQPointConversion(f, qpoint_face);

    const auto J = RefJacobian(qpoint);
    const auto JT = chi_math::Transpose(J);
    const auto JTinv = chi_math::Inverse(JT);
    const auto [detJ, qp_normal] =
      RefFaceJacobianDeterminantAndNormal(f, qpoint_face);

    normals[qp] = qp_normal;

    JxW[qp] = detJ * surface_quadrature.weights_[qp];

    for (size_t i = 0; i < num_nodes_; ++i)
    {
      const double ref_shape_i = RefShape(i, qpoint);
      shape_value[i][qp] = ref_shape_i;
      const Vec3 ref_shape_grad = RefGradShape(i, qpoint);
      const VecDbl b = {ref_shape_grad.x, ref_shape_grad.y, ref_shape_grad.z};
      const VecDbl JTInv_b = MatMul(JTinv, b);

      shape_grad[i][qp] = {JTInv_b[0], JTInv_b[1], JTInv_b[2]};

      const auto& node_xyz = node_locations_[i];

      qpoints_xyz[qp] += ref_shape_i * node_xyz;
    } // for i
  }   // for qp

  return finite_element::SurfaceQuadraturePointData(quadrature_point_indices,
                                                    qpoints_xyz,
                                                    shape_value,
                                                    shape_grad,
                                                    JxW,
                                                    normals,
                                                    face_node_mappings_,
                                                    num_nodes_);
}

std::pair<double, LagrangeBaseMapping::Vec3>
LagrangeBaseMapping::RefFaceJacobianDeterminantAndNormal(
  size_t face_index, const Vec3& qpoint_face) const
{
  ChiLogicalError("Method not implemented");
}

void LagrangeBaseMapping::ComputeCellVolumeAndAreas(
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
  else if (coordinate_system_type_ == CoordinateSystemType::SPHERICAL)
    swf = SpatialDiscretization::Spherical1DSpatialWeightFunction;
  else
    ChiInvalidArgument("Unsupported coordinate system encountered");

  volume = 0.0;
  {
    const size_t num_qpoints = volume_quadrature_.qpoints_.size();

    for (size_t qp = 0; qp < num_qpoints; ++qp)
    {
      const Vec3& qpoint = volume_quadrature_.qpoints_[qp];
      const auto J = RefJacobian(qpoint);
      const double detJ = Determinant(J);

      Vec3 qpoint_xyz = {0.0, 0.0, 0.0};
      for (size_t i = 0; i < num_nodes_; ++i)
        qpoint_xyz += RefShape(i, qpoint) * node_locations_[i];

      volume += swf(qpoint_xyz) * detJ * volume_quadrature_.weights_[qp];
    }
  } // volume

  areas.clear();
  {
    const size_t num_faces = cell.faces_.size();
    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face_quadrature = GetSurfaceQuadrature(f);
      const size_t num_qpoints = face_quadrature.qpoints_.size();

      double area = 0.0;
      for (size_t qp = 0; qp < num_qpoints; ++qp)
      {
        const Vec3& qpoint_face = face_quadrature.qpoints_[qp];
        const auto qpoint = FaceToElementQPointConversion(f, qpoint_face);
        const auto J = RefJacobian(qpoint);
        const auto [detJ, qp_normal] =
          RefFaceJacobianDeterminantAndNormal(f, qpoint_face);

        Vec3 qpoint_xyz = {0.0, 0.0, 0.0};
        for (size_t i = 0; i < num_nodes_; ++i)
          qpoint_xyz += RefShape(i, qpoint) * node_locations_[i];

        area += swf(qpoint_xyz) * detJ * face_quadrature.weights_[qp];
      }
      areas.push_back(area);
    } // for face
  }   // areas
}

// ##################################################################
WorldXYZToNaturalMappingHelper::WorldXYZToNaturalMappingHelper(
  const LagrangeBaseMapping& cell_mapping, const Vec3& world_x)
  : cell_mapping_(cell_mapping), world_x_(world_x)
{
  const auto& cell = cell_mapping.ReferenceCell();
  const auto cell_type = cell.Type();
  if (cell_type == chi_mesh::CellType::SLAB) dimension_ = 1;
  else if (cell_type == chi_mesh::CellType::POLYGON)
    dimension_ = 2;
  else if (cell_type == chi_mesh::CellType::POLYHEDRON)
    dimension_ = 3;
  else
    ChiLogicalError("Unsupported cell type.");
}

VecDbl WorldXYZToNaturalMappingHelper::F(const VecDbl& x) const
{
  const size_t num_nodes = cell_mapping_.NumNodes();
  const auto& node_x = cell_mapping_.GetNodeLocations();

  Vec3 result_vec3 = world_x_;
  if (dimension_ == 1)
    for (uint32_t i = 0; i < num_nodes; ++i)
    {
      const double shape_i_val = cell_mapping_.RefShape(i, {x[0], x[1], x[2]});
      result_vec3 -= Vec3(0.0, 0.0, shape_i_val * node_x[i].z);
    }
  else if (dimension_ == 2)
    for (uint32_t i = 0; i < num_nodes; ++i)
    {
      const double shape_i_val = cell_mapping_.RefShape(i, {x[0], x[1], x[2]});
      result_vec3 -=
        Vec3(shape_i_val * node_x[i].x, shape_i_val * node_x[i].y, 0.0);
    }
  else if (dimension_ == 3)
    for (uint32_t i = 0; i < num_nodes; ++i)
      result_vec3 -= cell_mapping_.RefShape(i, {x[0], x[1], x[2]}) * node_x[i];

  return {result_vec3.x, result_vec3.y, result_vec3.z};
}

MatDbl WorldXYZToNaturalMappingHelper::J(const VecDbl& x) const
{
  return cell_mapping_.RefJacobian({x[0], x[1], x[2]});
}

} // namespace chi_math::cell_mapping