#include "finite_element.h"

namespace chi_math
{
namespace finite_element
{
void InternalQuadraturePointData::InitializeData(
  std::vector<unsigned int> quadrature_point_indices,
  VecVec3 qpoints_xyz,
  std::vector<VecDbl> shape_value,
  std::vector<VecVec3> shape_grad,
  VecDbl JxW,
  std::vector<std::vector<int>> face_dof_mappings,
  size_t num_nodes)
{
  quadrature_point_indices_ = std::move(quadrature_point_indices);
  qpoints_xyz_ = std::move(qpoints_xyz);
  shape_value_ = std::move(shape_value);
  shape_grad_ = std::move(shape_grad);
  JxW_ = std::move(JxW);
  face_dof_mappings_ = std::move(face_dof_mappings);
  num_nodes_ = num_nodes;
  initialized_ = true;
}

void InternalQuadraturePointData::InitializeEmpty()
{
  initialized_ = true;
}

const std::vector<unsigned int>&
InternalQuadraturePointData::QuadraturePointIndices() const
{
  if (not initialized_) THROW_QP_UNINIT();
  return quadrature_point_indices_;
}
chi_mesh::Vector3 InternalQuadraturePointData::QPointXYZ(unsigned int qp) const
{
  if (not initialized_) THROW_QP_UNINIT();
  return qpoints_xyz_.at(qp);
}
double InternalQuadraturePointData::ShapeValue(unsigned int i,
                                               unsigned int qp) const
{
  if (not initialized_) THROW_QP_UNINIT();
  auto& qp_data = shape_value_.at(i);
  return qp_data.at(qp);
}
chi_mesh::Vector3 InternalQuadraturePointData::ShapeGrad(unsigned int i,
                                                         unsigned int qp) const
{
  if (not initialized_) THROW_QP_UNINIT();
  auto& qp_data = shape_grad_.at(i);
  return qp_data.at(qp);
}

const VecVec3& InternalQuadraturePointData::QPointsXYZ() const
{
  return qpoints_xyz_;
}

const std::vector<VecDbl>& InternalQuadraturePointData::ShapeValues() const
{
  return shape_value_;
}
const std::vector<VecVec3>& InternalQuadraturePointData::ShapeGradValues() const
{
  return shape_grad_;
}
const std::vector<double>& InternalQuadraturePointData::JxW_Values() const
{
  return JxW_;
}
double InternalQuadraturePointData::JxW(unsigned int qp) const
{
  if (not initialized_) THROW_QP_UNINIT();
  return JxW_.at(qp);
}
int InternalQuadraturePointData::FaceDofMapping(size_t face,
                                                size_t face_node_index) const
{
  if (not initialized_) THROW_QP_UNINIT();
  auto& face_data = face_dof_mappings_.at(face);
  return face_data.at(face_node_index);
}
size_t InternalQuadraturePointData::NumNodes() const
{
  if (not initialized_) THROW_QP_UNINIT();
  return num_nodes_;
}

void FaceQuadraturePointData::InitializeData(
  std::vector<unsigned int> quadrature_point_indices,
  VecVec3 qpoints_xyz,
  std::vector<VecDbl> shape_value,
  std::vector<VecVec3> shape_grad,
  VecDbl JxW,
  VecVec3 normals,
  std::vector<std::vector<int>> face_dof_mappings,
  size_t num_nodes)
{
  quadrature_point_indices_ = std::move(quadrature_point_indices);
  qpoints_xyz_ = std::move(qpoints_xyz);
  shape_value_ = std::move(shape_value);
  shape_grad_ = std::move(shape_grad);
  JxW_ = std::move(JxW);
  normals_ = std::move(normals);
  face_dof_mappings_ = std::move(face_dof_mappings);
  num_nodes_ = num_nodes;
  initialized_ = true;
}

chi_mesh::Vector3 FaceQuadraturePointData::Normal(unsigned int qp) const
{
  if (not initialized_) THROW_QP_UNINIT();
  return normals_.at(qp);
}

const VecVec3& FaceQuadraturePointData::Normals() const
{
  return normals_;
}
} // namespace finite_element
} // namespace chi_math
