#include "pwl_polygon.h"

void PolygonPWLFEValues::InitializeQuadraturePointData(
  chi_math::finite_element::InternalQuadraturePointData& internal_data,
  std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data)
{
  active_volume_quadrature = &arbitrary_volume_quadrature;
  active_surface_quadrature= &arbitrary_surface_quadrature;

  auto& vol_quadrature = arbitrary_volume_quadrature;
  auto& srf_quadrature = arbitrary_surface_quadrature;

  //=================================== Determine number of internal qpoints
  size_t num_tris = sides.size();
  size_t  num_vol_qpoints = vol_quadrature.qpoints.size();
  size_t ttl_num_vol_qpoints = num_tris * num_vol_qpoints;

  //=================================== Declare necessary vars
  std::vector<unsigned int>     V_quadrature_point_indices;
  VecVec3                       V_qpoints_xyz;
  std::vector<VecDbl>           V_shape_value;
  std::vector<VecVec3>          V_shape_grad;
  VecDbl                        V_JxW;
  size_t                        V_num_nodes;

  //=================================== Init volumetric quadrature
  V_quadrature_point_indices.reserve(ttl_num_vol_qpoints);
  V_qpoints_xyz.reserve(ttl_num_vol_qpoints);
  for (unsigned int qp=0; qp<ttl_num_vol_qpoints; ++qp)
    V_quadrature_point_indices.push_back(qp);

  V_shape_value.reserve(num_nodes);
  V_shape_grad.reserve(num_nodes);
  V_JxW.reserve(num_nodes);
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
        auto& v0= sides[s].v0;
        auto& J = sides[s].J;
        auto& qp_xyz_tilde = vol_quadrature.qpoints[qp];
        V_qpoints_xyz.push_back(v0 + J * qp_xyz_tilde);
      }//for qp
    } //for side

    V_shape_value.push_back(node_shape_value);
    V_shape_grad.push_back(node_shape_grad);
  }//for i

  V_JxW.reserve(ttl_num_vol_qpoints);
  for (auto& side : sides)
  {
    for (size_t qp=0; qp<num_vol_qpoints; ++qp)
    {
      double w = vol_quadrature.weights[qp];
      V_JxW.push_back(side.detJ * w);
    }//for qp
  } //for side
  V_num_nodes = num_nodes;
  internal_data.InitializeData(V_quadrature_point_indices,
                               V_qpoints_xyz,
                               V_shape_value,
                               V_shape_grad,
                               V_JxW,
                               face_dof_mappings,
                               V_num_nodes);

  //=================================== Init surface quadrature
  size_t num_srf_qpoints = srf_quadrature.qpoints.size();
  faces_qp_data.resize(sides.size());

  for (unsigned int s=0; s<sides.size(); ++s)
  {
    //=================================== Declare necessary vars
    std::vector<unsigned int>     F_quadrature_point_indices;
    VecVec3                       F_qpoints_xyz;
    std::vector<VecDbl>           F_shape_value;
    std::vector<VecVec3>          F_shape_grad;
    VecDbl                        F_JxW;
    VecVec3                       F_normals;
    size_t                        F_num_nodes;

    size_t ttl_num_face_qpoints = num_srf_qpoints;

    F_quadrature_point_indices.reserve(ttl_num_face_qpoints);
    for (unsigned int qp=0; qp<ttl_num_face_qpoints; ++qp)
      F_quadrature_point_indices.push_back(qp);

    for (size_t qp=0; qp<num_vol_qpoints; ++qp)
      F_normals.push_back(sides[s].normal);

    F_shape_value.reserve(num_nodes);
    F_shape_grad.reserve(num_nodes);
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
      F_shape_value.push_back(node_shape_value);
      F_shape_grad.push_back(node_shape_grad);
    }//for i

    F_JxW.reserve(ttl_num_face_qpoints);
    for (size_t qp=0; qp<num_srf_qpoints; ++qp)
    {
      double w = srf_quadrature.weights[qp];
      F_JxW.push_back(sides[s].detJ_surf * w);
    }
    F_num_nodes = 2;
    auto& face_qp_data = faces_qp_data[s];
    face_qp_data.InitializeData(F_quadrature_point_indices,
                                F_qpoints_xyz,
                                F_shape_value,
                                F_shape_grad,
                                F_JxW,
                                F_normals,
                                face_dof_mappings,
                                F_num_nodes);
  }//for face
}