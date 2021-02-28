#include "finite_element.h"

namespace chi_math
{
namespace finite_element
{
  void InternalQuadraturePointData::
       InitializeData(std::vector<unsigned int> quadrature_point_indices,
                      VecVec3                   qpoints_xyz,
                      std::vector<VecDbl>       shape_value,
                      std::vector<VecVec3>      shape_grad,
                      VecDbl                    JxW,
                      std::vector<std::vector<int>> face_dof_mappings,
                      size_t num_nodes)
  {
    m_quadrature_point_indices= std::move(quadrature_point_indices);
    m_qpoints_xyz             = std::move(qpoints_xyz             );
    m_shape_value             = std::move(shape_value             );
    m_shape_grad              = std::move(shape_grad              );
    m_JxW                     = std::move(JxW                     );
    m_face_dof_mappings       = std::move(face_dof_mappings       );
    m_num_nodes               = num_nodes               ;
    m_initialized = true;
  }

  const std::vector<unsigned int>& InternalQuadraturePointData::
    QuadraturePointIndices() const
  {
    if (not m_initialized) THROW_QP_UNINIT();
    return m_quadrature_point_indices;
  }
  chi_mesh::Vector3 InternalQuadraturePointData::
    QPointXYZ(unsigned int qp) const
  {
    if (not m_initialized) THROW_QP_UNINIT();
    return m_qpoints_xyz.at(qp);
  }
  double InternalQuadraturePointData::
    ShapeValue(unsigned int i, unsigned int qp) const
  {
    if (not m_initialized) THROW_QP_UNINIT();
    auto& qp_data = m_shape_value.at(i);
    return qp_data.at(qp);
  }
  chi_mesh::Vector3 InternalQuadraturePointData::
    ShapeGrad(unsigned int i, unsigned int qp) const
  {
    if (not m_initialized) THROW_QP_UNINIT();
    auto& qp_data = m_shape_grad.at(i);
    return qp_data.at(qp);
  }
  double InternalQuadraturePointData::JxW(unsigned int qp) const
  {
    if (not m_initialized) THROW_QP_UNINIT();
    return m_JxW.at(qp);
  }
  int InternalQuadraturePointData::FaceDofMapping(size_t face,
                                                  size_t face_node_index) const
  {
    if (not m_initialized) THROW_QP_UNINIT();
    auto& face_data = m_face_dof_mappings.at(face);
    return face_data.at(face_node_index);
  }
  size_t InternalQuadraturePointData::NumNodes() const
  {
    if (not m_initialized) THROW_QP_UNINIT();
    return m_num_nodes;
  }





  void FaceQuadraturePointData::
       InitializeData(std::vector<unsigned int>     quadrature_point_indices,
                      VecVec3                       qpoints_xyz,
                      std::vector<VecDbl>           shape_value,
                      std::vector<VecVec3>          shape_grad,
                      VecDbl                        JxW,
                      VecVec3                       normals,
                      std::vector<std::vector<int>> face_dof_mappings,
                      size_t num_nodes)
  {
    m_quadrature_point_indices= std::move(quadrature_point_indices);
    m_qpoints_xyz             = std::move(qpoints_xyz             );
    m_shape_value             = std::move(shape_value             );
    m_shape_grad              = std::move(shape_grad              );
    m_JxW                     = std::move(JxW                     );
    m_normals                 = std::move(normals                 );
    m_face_dof_mappings       = std::move(face_dof_mappings       );
    m_num_nodes               = num_nodes               ;
    m_initialized = true;
  }

  double FaceQuadraturePointData::Normal(unsigned int qp) const
  {
    if (not m_initialized) THROW_QP_UNINIT();
    return m_JxW.at(qp);
  }
}
}

