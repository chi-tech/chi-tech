#include "pwl_polyhedron.h"

//###################################################################
/**Computes cell volume and surface integrals.*/
void PolyhedronPWLFEValues::ComputeUnitIntegrals()
{
  const bool ON_SURFACE = true;
  if (precomputed) return;

  //============================================= Precompute elements
  for (size_t f=0; f < face_data.size(); f++)
  {
    for (size_t s=0; s < face_data[f].sides.size(); s++)
    {
      for (size_t i=0; i < num_nodes; i++)
      {
        FEqp_data3d pernode_data;

        size_t num_vol_qpoints = default_volume_quadrature.qpoints.size();
        size_t num_srf_qpoints = default_surface_quadrature.qpoints.size();

        //========== Reserving data
        pernode_data.gradshapex_qp.reserve(num_vol_qpoints);
        pernode_data.gradshapey_qp.reserve(num_vol_qpoints);
        pernode_data.gradshapez_qp.reserve(num_vol_qpoints);
        pernode_data.shape_qp.reserve(num_vol_qpoints);
        pernode_data.shape_qp_surf.reserve(num_srf_qpoints);

        //Prestore GradVarphi_xyz
        for (size_t qp=0; qp<num_vol_qpoints; qp++)
        {
          pernode_data.gradshapex_qp.push_back(FaceSideGradShape_x(f, s, i));
          pernode_data.gradshapey_qp.push_back(FaceSideGradShape_y(f, s, i));
          pernode_data.gradshapez_qp.push_back(FaceSideGradShape_z(f, s, i));
        }
        //Prestore Varphi
        for (size_t qp=0; qp<num_vol_qpoints; qp++)
          pernode_data.shape_qp.push_back(FaceSideShape(f, s, i, qp));

        //Prestore Varphi on surface
        for (size_t qp=0; qp<num_srf_qpoints; qp++)
          pernode_data.shape_qp_surf.push_back(FaceSideShape(f, s, i, qp, ON_SURFACE));

        face_data[f].sides[s].qp_data.push_back(pernode_data);
      } // for i

    } //for side
  } //for face

  //============================================= Lambdas for accessing data
  /**Determinant evaluated at quadrature point*/
  auto DetJ = [this](int face_index, int side_index,
                     int qpoint_index, bool on_surface=false)
  {
    if (on_surface)
      return (face_data[face_index].sides[side_index].detJ_surf);
    else
      return (face_data[face_index].sides[side_index].detJ);
  };

  /**Shape function evaluation on a tet at a quadrature point*/
  auto GetShape = [this](int face, int side, int i, int qp, bool surface = false)
  {
    if (surface)
      return face_data[face].sides[side].qp_data[i].shape_qp_surf[qp];
    else
      return face_data[face].sides[side].qp_data[i].shape_qp[qp];
  };

  /**GradeShape-x function evaluation on a tet at a quadrature point*/
  auto GetGradShape_x = [this](int face, int side, int i, int qp)
  { return face_data[face].sides[side].qp_data[i].gradshapex_qp[qp]; };

  /**GradeShape-y function evaluation on a tet at a quadrature point*/
  auto GetGradShape_y = [this](int face, int side, int i, int qp)
  { return face_data[face].sides[side].qp_data[i].gradshapey_qp[qp]; };

  /**GradeShape-z function evaluation on a tet at a quadrature point*/
  auto GetGradShape_z = [this](int face, int side, int i, int qp)
  { return face_data[face].sides[side].qp_data[i].gradshapez_qp[qp]; };

  //============================================= Volume integrals
  IntV_gradShapeI_gradShapeJ.reserve(num_nodes);
  IntV_shapeI_gradshapeJ.reserve(num_nodes);
  IntV_shapeI_shapeJ.reserve(num_nodes);
  IntV_shapeI.reserve(num_nodes);

  for (int i=0; i < num_nodes; i++)
  {
    std::vector<double> gradijvalue_i(num_nodes, 0.0);
    std::vector<chi_mesh::Vector3> varphi_i_gradj(num_nodes, chi_mesh::Vector3());
    std::vector<double>           varphi_i_varphi_j(num_nodes, 0);

    // Computing
    // GradVarphi_i*GradVarphi_j and
    // Varphi_i*GradVarphi_j_xyz and
    // Varphi_i*Varphi_j
    for (int j=0; j < num_nodes; j++)
    {
      gradijvalue_i[j] = 0.0;
      varphi_i_varphi_j[j] = 0.0;

      for (size_t f=0; f < face_data.size(); f++)
      {
        for (size_t s = 0; s < face_data[f].sides.size(); s++)
        {
          for (size_t qp=0; qp < default_volume_quadrature.qpoints.size(); qp++)
          {
            gradijvalue_i[j]
              += default_volume_quadrature.weights[qp] *
                 (GetGradShape_x(f, s, i, qp)*GetGradShape_x(f, s, j, qp) +
                  GetGradShape_y(f, s, i, qp)*GetGradShape_y(f, s, j, qp) +
                  GetGradShape_z(f, s, i, qp)*GetGradShape_z(f, s, j, qp)) *
                 DetJ(f,s,qp);
          }//for qp

          for (size_t qp=0; qp < default_volume_quadrature.qpoints.size(); qp++)
          {
            double varphi_i = GetShape(f, s, i, qp);
            double weight = default_volume_quadrature.weights[qp];

            varphi_i_gradj[j].x
              += weight*varphi_i* GetGradShape_x(f, s, j, qp)*DetJ(f,s,qp);

            varphi_i_gradj[j].y
              += weight*varphi_i* GetGradShape_y(f, s, j, qp)*DetJ(f,s,qp);

            varphi_i_gradj[j].z
              += weight*varphi_i* GetGradShape_z(f, s, j, qp)*DetJ(f,s,qp);

            varphi_i_varphi_j[j]
              += weight*varphi_i* GetShape(f, s, j, qp)*DetJ(f,s,qp);
          }// for qp
        }// for s
      }// for f
    }// for j

    //Computing Varphi_i
    double  valuei_i = 0.0;
    for (size_t f=0; f < face_data.size(); f++)
    {
      for (size_t s = 0; s < face_data[f].sides.size(); s++)
      {
        for (size_t qp=0; qp < default_volume_quadrature.qpoints.size(); qp++)
        {
          valuei_i += default_volume_quadrature.weights[qp] *
                      GetShape(f, s, i, qp) *
                      DetJ(f,s,qp);
        }// for gp
      } // for s
    }// for f
    IntV_gradShapeI_gradShapeJ.push_back(std::move(gradijvalue_i));
    IntV_shapeI_gradshapeJ.push_back(std::move(varphi_i_gradj));
    IntV_shapeI_shapeJ.push_back(std::move(varphi_i_varphi_j));
    IntV_shapeI.push_back(valuei_i);
  }

  //============================================= Surface integrals
  std::vector<std::vector<std::vector<double>>> IntSi_shapeI_shapeJ;
  std::vector<std::vector<std::vector<chi_mesh::Vector3>>> IntSi_shapeI_gradshapeJ;

  for (int i=0; i < num_nodes; i++)
  {
    // Computing
    // Varphi_i*Varphi_j on each face and
    // Varphi_i on each face
    std::vector<std::vector<double>> varphi_i_varphi_j_surf;
    std::vector<std::vector<chi_mesh::Vector3>> varphi_i_gradvarphi_j_surf;
    std::vector<double> varphi_i_surf(face_data.size(), 0.0);

    for (size_t f=0; f < face_data.size(); f++)
    {
      std::vector<double> f_varphi_i_varphi_j_surf(num_nodes, 0);
      std::vector<chi_mesh::Vector3> f_varphi_i_grad_j_surf;
      f_varphi_i_grad_j_surf.resize(num_nodes, chi_mesh::Vector3());

      for (int j=0; j < num_nodes; j++)
      {
        double value_ij = 0.0;
        double value_x_ij = 0.0;
        double value_y_ij = 0.0;
        double value_z_ij = 0.0;

        for (size_t s = 0; s < face_data[f].sides.size(); s++)
        {
          for (size_t qp=0; qp < default_surface_quadrature.qpoints.size(); qp++)
          {
            value_ij
              += default_surface_quadrature.weights[qp] *
                 GetShape(f, s, i, qp, ON_SURFACE) *
                 GetShape(f, s, j, qp, ON_SURFACE) *
                 DetJ(f,s,qp,ON_SURFACE);

            value_x_ij
              += default_surface_quadrature.weights[qp] *
                 GetShape(f, s, i, qp, ON_SURFACE) *
                 GetGradShape_x(f, s, j, qp) *
                 DetJ(f,s,qp,ON_SURFACE);

            value_y_ij
              += default_surface_quadrature.weights[qp] *
                 GetShape(f, s, i, qp, ON_SURFACE) *
                 GetGradShape_y(f, s, j, qp) *
                 DetJ(f,s,qp,ON_SURFACE);

            value_z_ij
              += default_surface_quadrature.weights[qp] *
                 GetShape(f, s, i, qp, ON_SURFACE) *
                 GetGradShape_z(f, s, j, qp) *
                 DetJ(f,s,qp,ON_SURFACE);
          }// for gp
        }

        f_varphi_i_varphi_j_surf[j] = value_ij;
        f_varphi_i_grad_j_surf[j].x = value_x_ij;
        f_varphi_i_grad_j_surf[j].y = value_y_ij;
        f_varphi_i_grad_j_surf[j].z = value_z_ij;
      }// for j
      varphi_i_varphi_j_surf.push_back(f_varphi_i_varphi_j_surf);
      varphi_i_gradvarphi_j_surf.push_back(f_varphi_i_grad_j_surf);

      double f_varphi_i_surf = 0.0;
      for (size_t s = 0; s < face_data[f].sides.size(); s++)
      {
        for (size_t qp=0; qp < default_surface_quadrature.qpoints.size(); qp++)
        {
          f_varphi_i_surf
            += default_surface_quadrature.weights[qp] *
               GetShape(f, s, i, qp, ON_SURFACE) *
               DetJ(f,s,qp,ON_SURFACE);
        }// for gp
      }
      varphi_i_surf[f] = f_varphi_i_surf;
    }// for f
    IntSi_shapeI_shapeJ.push_back(std::move(varphi_i_varphi_j_surf));
    IntSi_shapeI_gradshapeJ.push_back(std::move(varphi_i_gradvarphi_j_surf));
    IntS_shapeI.push_back(std::move(varphi_i_surf));
  }// for i

  //============================================= Reindexing surface integrals
  IntS_shapeI_shapeJ.resize(face_data.size());
  IntS_shapeI_gradshapeJ.resize(face_data.size());
  for (size_t f=0; f < face_data.size(); f++)
  {
    IntS_shapeI_shapeJ[f].resize(num_nodes);
    IntS_shapeI_gradshapeJ[f].resize(num_nodes);
    for (int i=0; i < num_nodes; i++)
    {
      IntS_shapeI_shapeJ[f][i].resize(num_nodes);
      IntS_shapeI_gradshapeJ[f][i].resize(num_nodes);
      for (int j=0; j < num_nodes; j++)
      {
        IntS_shapeI_shapeJ[f][i][j] = IntSi_shapeI_shapeJ[i][f][j];
        IntS_shapeI_gradshapeJ[f][i][j] = IntSi_shapeI_gradshapeJ[i][f][j];
      }
    }
  }

  //============================================= Cleanup
  for (auto& face : face_data)
    for (auto& side : face.sides)
      side.qp_data = std::vector<FEqp_data3d>(0);

  precomputed = true;
}

//###################################################################
/**Computes cell volume and surface integrals.*/
void PolyhedronPWLFEValues::
  ComputeUnitIntegrals(chi_math::finite_element::UnitIntegralData& ui_data)
{
  const bool ON_SURFACE = true;
  if (precomputed) return;

  //============================================= Precompute elements
  for (size_t f=0; f < face_data.size(); f++)
  {
    for (size_t s=0; s < face_data[f].sides.size(); s++)
    {
      for (size_t i=0; i < num_nodes; i++)
      {
        FEqp_data3d pernode_data;

        size_t num_vol_qpoints = default_volume_quadrature.qpoints.size();
        size_t num_srf_qpoints = default_surface_quadrature.qpoints.size();

        //========== Reserving data
        pernode_data.gradshapex_qp.reserve(num_vol_qpoints);
        pernode_data.gradshapey_qp.reserve(num_vol_qpoints);
        pernode_data.gradshapez_qp.reserve(num_vol_qpoints);
        pernode_data.shape_qp.reserve(num_vol_qpoints);
        pernode_data.shape_qp_surf.reserve(num_srf_qpoints);

        //Prestore GradVarphi_xyz
        for (size_t qp=0; qp<num_vol_qpoints; qp++)
        {
          pernode_data.gradshapex_qp.push_back(FaceSideGradShape_x(f, s, i));
          pernode_data.gradshapey_qp.push_back(FaceSideGradShape_y(f, s, i));
          pernode_data.gradshapez_qp.push_back(FaceSideGradShape_z(f, s, i));
        }
        //Prestore Varphi
        for (size_t qp=0; qp<num_vol_qpoints; qp++)
          pernode_data.shape_qp.push_back(FaceSideShape(f, s, i, qp));

        //Prestore Varphi on surface
        for (size_t qp=0; qp<num_srf_qpoints; qp++)
          pernode_data.shape_qp_surf.push_back(FaceSideShape(f, s, i, qp, ON_SURFACE));

        face_data[f].sides[s].qp_data.push_back(pernode_data);
      } // for i

    } //for side
  } //for face

  //============================================= Lambdas for accessing data
  /**Determinant evaluated at quadrature point*/
  auto DetJ = [this](int face_index, int side_index,
                     int qpoint_index, bool on_surface=false)
  {
    if (on_surface)
      return (face_data[face_index].sides[side_index].detJ_surf);
    else
      return (face_data[face_index].sides[side_index].detJ);
  };

  /**Shape function evaluation on a tet at a quadrature point*/
  auto GetShape = [this](int face, int side, int i, int qp, bool surface = false)
  {
    if (surface)
      return face_data[face].sides[side].qp_data[i].shape_qp_surf[qp];
    else
      return face_data[face].sides[side].qp_data[i].shape_qp[qp];
  };

  /**GradeShape-x function evaluation on a tet at a quadrature point*/
  auto GetGradShape_x = [this](int face, int side, int i, int qp)
  { return face_data[face].sides[side].qp_data[i].gradshapex_qp[qp]; };

  /**GradeShape-y function evaluation on a tet at a quadrature point*/
  auto GetGradShape_y = [this](int face, int side, int i, int qp)
  { return face_data[face].sides[side].qp_data[i].gradshapey_qp[qp]; };

  /**GradeShape-z function evaluation on a tet at a quadrature point*/
  auto GetGradShape_z = [this](int face, int side, int i, int qp)
  { return face_data[face].sides[side].qp_data[i].gradshapez_qp[qp]; };

  //============================================= Volume integrals
  ui_data.IntV_gradShapeI_gradShapeJ.reserve(num_nodes);
  ui_data.IntV_shapeI_gradshapeJ.reserve(num_nodes);
  ui_data.IntV_shapeI_shapeJ.reserve(num_nodes);
  ui_data.IntV_shapeI.reserve(num_nodes);

  for (int i=0; i < num_nodes; i++)
  {
    std::vector<double> gradijvalue_i(num_nodes, 0.0);
    std::vector<chi_mesh::Vector3> varphi_i_gradj(num_nodes, chi_mesh::Vector3());
    std::vector<double>           varphi_i_varphi_j(num_nodes, 0);

    // Computing
    // GradVarphi_i*GradVarphi_j and
    // Varphi_i*GradVarphi_j_xyz and
    // Varphi_i*Varphi_j
    for (int j=0; j < num_nodes; j++)
    {
      gradijvalue_i[j] = 0.0;
      varphi_i_varphi_j[j] = 0.0;

      for (size_t f=0; f < face_data.size(); f++)
      {
        for (size_t s = 0; s < face_data[f].sides.size(); s++)
        {
          for (size_t qp=0; qp < default_volume_quadrature.qpoints.size(); qp++)
          {
            gradijvalue_i[j]
              += default_volume_quadrature.weights[qp] *
                 (GetGradShape_x(f, s, i, qp)*GetGradShape_x(f, s, j, qp) +
                  GetGradShape_y(f, s, i, qp)*GetGradShape_y(f, s, j, qp) +
                  GetGradShape_z(f, s, i, qp)*GetGradShape_z(f, s, j, qp)) *
                 DetJ(f,s,qp);
          }//for qp

          for (size_t qp=0; qp < default_volume_quadrature.qpoints.size(); qp++)
          {
            double varphi_i = GetShape(f, s, i, qp);
            double weight = default_volume_quadrature.weights[qp];

            varphi_i_gradj[j].x
              += weight*varphi_i* GetGradShape_x(f, s, j, qp)*DetJ(f,s,qp);

            varphi_i_gradj[j].y
              += weight*varphi_i* GetGradShape_y(f, s, j, qp)*DetJ(f,s,qp);

            varphi_i_gradj[j].z
              += weight*varphi_i* GetGradShape_z(f, s, j, qp)*DetJ(f,s,qp);

            varphi_i_varphi_j[j]
              += weight*varphi_i* GetShape(f, s, j, qp)*DetJ(f,s,qp);
          }// for qp
        }// for s
      }// for f
    }// for j

    //Computing Varphi_i
    double  valuei_i = 0.0;
    for (size_t f=0; f < face_data.size(); f++)
    {
      for (size_t s = 0; s < face_data[f].sides.size(); s++)
      {
        for (size_t qp=0; qp < default_volume_quadrature.qpoints.size(); qp++)
        {
          valuei_i += default_volume_quadrature.weights[qp] *
                      GetShape(f, s, i, qp) *
                      DetJ(f,s,qp);
        }// for gp
      } // for s
    }// for f
    ui_data.IntV_gradShapeI_gradShapeJ.push_back(std::move(gradijvalue_i));
    ui_data.IntV_shapeI_gradshapeJ.push_back(std::move(varphi_i_gradj));
    ui_data.IntV_shapeI_shapeJ.push_back(std::move(varphi_i_varphi_j));
    ui_data.IntV_shapeI.push_back(valuei_i);
  }

  //============================================= Surface integrals
  std::vector<std::vector<std::vector<double>>> IntSi_shapeI_shapeJ;
  std::vector<std::vector<std::vector<chi_mesh::Vector3>>> IntSi_shapeI_gradshapeJ;

  for (int i=0; i < num_nodes; i++)
  {
    // Computing
    // Varphi_i*Varphi_j on each face and
    // Varphi_i on each face
    std::vector<std::vector<double>> varphi_i_varphi_j_surf;
    std::vector<std::vector<chi_mesh::Vector3>> varphi_i_gradvarphi_j_surf;
    std::vector<double> varphi_i_surf(face_data.size(), 0.0);

    for (size_t f=0; f < face_data.size(); f++)
    {
      std::vector<double> f_varphi_i_varphi_j_surf(num_nodes, 0);
      std::vector<chi_mesh::Vector3> f_varphi_i_grad_j_surf;
      f_varphi_i_grad_j_surf.resize(num_nodes, chi_mesh::Vector3());

      for (int j=0; j < num_nodes; j++)
      {
        double value_ij = 0.0;
        double value_x_ij = 0.0;
        double value_y_ij = 0.0;
        double value_z_ij = 0.0;

        for (size_t s = 0; s < face_data[f].sides.size(); s++)
        {
          for (size_t qp=0; qp < default_surface_quadrature.qpoints.size(); qp++)
          {
            value_ij
              += default_surface_quadrature.weights[qp] *
                 GetShape(f, s, i, qp, ON_SURFACE) *
                 GetShape(f, s, j, qp, ON_SURFACE) *
                 DetJ(f,s,qp,ON_SURFACE);

            value_x_ij
              += default_surface_quadrature.weights[qp] *
                 GetShape(f, s, i, qp, ON_SURFACE) *
                 GetGradShape_x(f, s, j, qp) *
                 DetJ(f,s,qp,ON_SURFACE);

            value_y_ij
              += default_surface_quadrature.weights[qp] *
                 GetShape(f, s, i, qp, ON_SURFACE) *
                 GetGradShape_y(f, s, j, qp) *
                 DetJ(f,s,qp,ON_SURFACE);

            value_z_ij
              += default_surface_quadrature.weights[qp] *
                 GetShape(f, s, i, qp, ON_SURFACE) *
                 GetGradShape_z(f, s, j, qp) *
                 DetJ(f,s,qp,ON_SURFACE);
          }// for gp
        }

        f_varphi_i_varphi_j_surf[j] = value_ij;
        f_varphi_i_grad_j_surf[j].x = value_x_ij;
        f_varphi_i_grad_j_surf[j].y = value_y_ij;
        f_varphi_i_grad_j_surf[j].z = value_z_ij;
      }// for j
      varphi_i_varphi_j_surf.push_back(f_varphi_i_varphi_j_surf);
      varphi_i_gradvarphi_j_surf.push_back(f_varphi_i_grad_j_surf);

      double f_varphi_i_surf = 0.0;
      for (size_t s = 0; s < face_data[f].sides.size(); s++)
      {
        for (size_t qp=0; qp < default_surface_quadrature.qpoints.size(); qp++)
        {
          f_varphi_i_surf
            += default_surface_quadrature.weights[qp] *
               GetShape(f, s, i, qp, ON_SURFACE) *
               DetJ(f,s,qp,ON_SURFACE);
        }// for gp
      }
      varphi_i_surf[f] = f_varphi_i_surf;
    }// for f
    IntSi_shapeI_shapeJ.push_back(std::move(varphi_i_varphi_j_surf));
    IntSi_shapeI_gradshapeJ.push_back(std::move(varphi_i_gradvarphi_j_surf));
    ui_data.IntS_shapeI.push_back(std::move(varphi_i_surf));
  }// for i

  //============================================= Reindexing surface integrals
  ui_data.IntS_shapeI_shapeJ.resize(face_data.size());
  ui_data.IntS_shapeI_gradshapeJ.resize(face_data.size());
  for (size_t f=0; f < face_data.size(); f++)
  {
    ui_data.IntS_shapeI_shapeJ[f].resize(num_nodes);
    ui_data.IntS_shapeI_gradshapeJ[f].resize(num_nodes);
    for (int i=0; i < num_nodes; i++)
    {
      ui_data.IntS_shapeI_shapeJ[f][i].resize(num_nodes);
      ui_data.IntS_shapeI_gradshapeJ[f][i].resize(num_nodes);
      for (int j=0; j < num_nodes; j++)
      {
        ui_data.IntS_shapeI_shapeJ[f][i][j] = IntSi_shapeI_shapeJ[i][f][j];
        ui_data.IntS_shapeI_gradshapeJ[f][i][j] = IntSi_shapeI_gradshapeJ[i][f][j];
      }
    }
  }

  //============================================= Cleanup
  for (auto& face : face_data)
    for (auto& side : face.sides)
      side.qp_data = std::vector<FEqp_data3d>(0);

  precomputed = true;
}

/**Precomputes cell volume and surface integrals.*/
void PolyhedronPWLFEValues::PreComputeValues()
{
  ComputeUnitIntegrals();
}