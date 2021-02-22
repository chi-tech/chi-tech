#include "pwl_polygon.h"

void PolygonPWLFEValues::InitializeQuadraturePointData(
  chi_math::finite_element::InternalQuadraturePointData& internal_data,
  std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data)
{
  auto& vol_quadrature = default_volume_quadrature;
  auto& srf_quadrature = default_surface_quadrature;

  //=================================== Determine number of internal qpoints
  size_t num_tris = sides.size();
  size_t  num_vol_qpoints = vol_quadrature.qpoints.size();
  size_t ttl_num_vol_qpoints = num_tris * num_vol_qpoints;

  //=================================== Init volumetric quadrature
  internal_data.quadrature_point_indices.reserve(ttl_num_vol_qpoints);
  for (unsigned int qp=0; qp<ttl_num_vol_qpoints; ++qp)
    internal_data.quadrature_point_indices.push_back(qp);

  internal_data.m_shape_value.reserve(num_nodes);
  internal_data.m_shape_grad.reserve(num_nodes);
  internal_data.m_JxW.reserve(num_nodes);
  for (size_t i=0; i < num_nodes; i++)
  {
    VecDbl  node_shape_value;
    VecVec3 node_shape_grad;

    node_shape_value.reserve(ttl_num_vol_qpoints);
    node_shape_grad.reserve(ttl_num_vol_qpoints);

    for (size_t s=0; s < sides.size(); s++)
    {
      for (size_t qp=0; qp<num_vol_qpoints; ++qp)
      {
        node_shape_value.push_back(SideShape(s,i,qp));
        node_shape_grad.emplace_back(SideGradShape_x(s,i), //x
                                     SideGradShape_y(s,i), //y
                                     0.0);                 //z
      }//for qp
    } //for side

    internal_data.m_shape_value.push_back(node_shape_value);
    internal_data.m_shape_grad.push_back(node_shape_grad);
  }//for i

  internal_data.m_JxW.reserve(ttl_num_vol_qpoints);
  for (auto& side : sides)
  {
    for (size_t qp=0; qp<num_vol_qpoints; ++qp)
    {
      double w = vol_quadrature.weights[qp];
      internal_data.m_JxW.push_back(side.detJ * w);
    }//for qp
  } //for side
  internal_data.face_dof_mappings = face_dof_mappings;
  internal_data.num_nodes = num_nodes;
  internal_data.initialized = true;

  //=================================== Init surface quadrature
  size_t num_srf_qpoints = srf_quadrature.qpoints.size();
  faces_qp_data.resize(sides.size());

  for (unsigned int s=0; s<sides.size(); ++s)
  {
    auto& face_qp_data = faces_qp_data[s];

    size_t ttl_num_face_qpoints = num_srf_qpoints;

    face_qp_data.quadrature_point_indices.reserve(ttl_num_face_qpoints);
    for (unsigned int qp=0; qp<ttl_num_face_qpoints; ++qp)
      face_qp_data.quadrature_point_indices.push_back(qp);

    for (size_t qp=0; qp<num_vol_qpoints; ++qp)
      face_qp_data.m_normals.push_back(sides[s].normal);

    face_qp_data.m_shape_value.reserve(num_nodes);
    face_qp_data.m_shape_grad.reserve(num_nodes);
    for (size_t i=0; i < num_nodes; i++)
    {
      VecDbl  node_shape_value;
      VecVec3 node_shape_grad;

      node_shape_value.reserve(ttl_num_face_qpoints);
      node_shape_grad.reserve(ttl_num_face_qpoints);

      for (size_t qp=0; qp<num_srf_qpoints; ++qp)
      {
        node_shape_value.push_back(SideShape(s,i,qp,true));
        node_shape_grad.emplace_back(SideGradShape_x(s,i), //x
                                     SideGradShape_y(s,i), //y
                                     0.0);                 //z
      }//for qp
      face_qp_data.m_shape_value.push_back(node_shape_value);
      face_qp_data.m_shape_grad.push_back(node_shape_grad);
    }//for i

    face_qp_data.m_JxW.reserve(ttl_num_face_qpoints);
    for (size_t qp=0; qp<num_srf_qpoints; ++qp)
    {
      double w = srf_quadrature.weights[qp];
      face_qp_data.m_JxW.push_back(sides[s].detJ_surf * w);
    }
    face_qp_data.face_dof_mappings = face_dof_mappings;
    face_qp_data.num_nodes = 2;
    face_qp_data.initialized = true;
  }//for face
}