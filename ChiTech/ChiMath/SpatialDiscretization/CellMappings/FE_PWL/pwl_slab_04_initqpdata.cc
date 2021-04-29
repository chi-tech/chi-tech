#include "pwl_slab.h"

void SlabMappingFE_PWL::InitializeAllQuadraturePointData(
  chi_math::finite_element::InternalQuadraturePointData& internal_data,
  std::vector<chi_math::finite_element::FaceQuadraturePointData>& faces_qp_data)
{
  InitializeVolumeQuadraturePointData(internal_data);
  faces_qp_data.resize(2);
  InitializeVolumeQuadraturePointData(faces_qp_data[0]);
  InitializeVolumeQuadraturePointData(faces_qp_data[1]);
}

void SlabMappingFE_PWL::InitializeVolumeQuadraturePointData(
  chi_math::finite_element::InternalQuadraturePointData& internal_data)
{
  auto& volume_quadrature = arbitrary_volume_quadrature;

  //=================================== Determine number of internal qpoints
  size_t ttl_num_vol_qpoints = volume_quadrature.qpoints.size();

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

    for (size_t qp=0; qp<ttl_num_vol_qpoints; ++qp)
    {
      node_shape_value.push_back(SlabShape(i,qp));
      node_shape_grad.emplace_back(0.0,                //x
                                   0.0,                //y
                                   SlabGradShape(i));  //z

      const double qp_xyz_tilde = volume_quadrature.qpoints[qp][0];
      V_qpoints_xyz.push_back(v0 +
                              h*chi_mesh::Vector3(0.0,0.0,qp_xyz_tilde));
    }//for qp

    V_shape_value.push_back(node_shape_value);
    V_shape_grad.push_back(node_shape_grad);
  }//for i

  V_JxW.reserve(ttl_num_vol_qpoints);
  const double J = h;
  for (size_t qp=0; qp<ttl_num_vol_qpoints; ++qp)
  {
    const double w = volume_quadrature.weights[qp];
    V_JxW.push_back(J * w);
  }//for qp
  V_num_nodes = num_nodes;
  internal_data.InitializeData(std::move(V_quadrature_point_indices),
                               std::move(V_qpoints_xyz             ),
                               std::move(V_shape_value             ),
                               std::move(V_shape_grad              ),
                               std::move(V_JxW                     ),
                               face_dof_mappings,
                               V_num_nodes);
}

void SlabMappingFE_PWL::InitializeFaceQuadraturePointData(unsigned int face,
                                                          chi_math::finite_element::FaceQuadraturePointData& faces_qp_data)
{

  //=================================== Init surface quadrature
  size_t num_srf_qpoints = 1;

  unsigned int f=face;
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

    for (size_t qp=0; qp<ttl_num_face_qpoints; ++qp)
      F_normals.push_back(normals[f]);

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
        node_shape_value.push_back(SlabShape(i,qp,true));
        node_shape_grad.emplace_back(0.0,                //x
                                     0.0,                //y
                                     SlabGradShape(i));  //z
      }//for qp
      F_shape_value.push_back(node_shape_value);
      F_shape_grad.push_back(node_shape_grad);
    }//for i

    F_JxW.reserve(ttl_num_face_qpoints);
    for (size_t qp=0; qp<num_srf_qpoints; ++qp)
    {
      double w = 1.0;
      F_JxW.push_back(1.0 * w);
    }
    F_num_nodes = 1;
    faces_qp_data.InitializeData(std::move(F_quadrature_point_indices),
                                 std::move(F_qpoints_xyz             ),
                                 std::move(F_shape_value             ),
                                 std::move(F_shape_grad              ),
                                 std::move(F_JxW                     ),
                                 std::move(F_normals                 ),
                                 face_dof_mappings,
                                 F_num_nodes);

  }//for face
}
