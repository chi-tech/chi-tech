#include "pwl_slab.h"

void SlabPWLFEView::InitializeQuadraturePointData()
{
  auto& vol_quadrature = default_volume_quadrature;

  //=================================== Determine number of internal qpoints
  size_t ttl_num_vol_qpoints = vol_quadrature.abscissae.size();

  //=================================== Init volumetric quadrature
  internal_qp_data.quadrature_point_indices.reserve(ttl_num_vol_qpoints);
  for (unsigned int qp=0; qp<ttl_num_vol_qpoints; ++qp)
    internal_qp_data.quadrature_point_indices.push_back(qp);

  internal_qp_data.m_shape_value.reserve(dofs);
  internal_qp_data.m_shape_grad.reserve(dofs);
  internal_qp_data.m_JxW.reserve(dofs);
  for (size_t i=0; i<dofs; i++)
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
    }//for qp

    internal_qp_data.m_shape_value.push_back(node_shape_value);
    internal_qp_data.m_shape_grad.push_back(node_shape_grad);
  }//for i

  internal_qp_data.m_JxW.reserve(ttl_num_vol_qpoints);
  double J = h/2.0;
  for (size_t qp=0; qp<ttl_num_vol_qpoints; ++qp)
  {
    double w = vol_quadrature.weights[qp];
    internal_qp_data.m_JxW.push_back(J * w);
  }//for qp
  internal_qp_data.initialized = true;

  //=================================== Init surface quadrature
  size_t num_srf_qpoints = 1;
  surface_qp_data.resize(2);

  for (unsigned int f=0; f<2; ++f)
  {
    auto& face_qp_data = surface_qp_data[f];

    size_t ttl_num_face_qpoints = num_srf_qpoints;

    face_qp_data.quadrature_point_indices.reserve(ttl_num_face_qpoints);
    for (unsigned int qp=0; qp<ttl_num_face_qpoints; ++qp)
      face_qp_data.quadrature_point_indices.push_back(qp);

    for (size_t qp=0; qp<ttl_num_face_qpoints; ++qp)
      face_qp_data.m_normals.push_back(normals[f]);

    face_qp_data.m_shape_value.reserve(dofs);
    face_qp_data.m_shape_grad.reserve(dofs);
    for (size_t i=0; i<dofs; i++)
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
      face_qp_data.m_shape_value.push_back(node_shape_value);
      face_qp_data.m_shape_grad.push_back(node_shape_grad);
    }//for i

    face_qp_data.m_JxW.reserve(ttl_num_face_qpoints);
    for (size_t qp=0; qp<num_srf_qpoints; ++qp)
    {
      double w = 1.0;
      face_qp_data.m_JxW.push_back(1.0 * w);
    }

    face_qp_data.initialized = true;
  }//for face
}