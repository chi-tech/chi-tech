#include "pwl_polygon.h"

#define ON_SURFACE true
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Varphi_x
/**Precomputation of the shape function at a quadrature point.*/
double PolygonPWLFEValues::PreShape(int s, int i, int qpoint_index, bool on_surface)
{
  double xi  = 0.0;
  double eta = 0.0;
  if (!on_surface)
  {
    auto& qpoint = volume_quadrature.qpoints.at(qpoint_index);

    xi = qpoint.x;
    eta= qpoint.y;
  }
  else
  {
    xi = surface_quadrature.qpoints[qpoint_index].x;
    eta = 0.0;
  }


  int index = node_to_side_map[i][s];
  double value = 0;
  if (index==0)
  {
    value = 1.0 - xi - eta;
  }
  if (index==1)
  {
    value = xi;
  }

  value += beta*eta;

  return value;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GradVarphi_x
/**Precomputation of the partial derivative along x of the
 * shape function at a quadrature point.*/
double PolygonPWLFEValues::PreGradShape_x(int s, int i, int qpoint_index)
{
  int index = node_to_side_map[i][s];
  double value = 0;
  if (index==0)
  {

    value = sides[s]->JTinv.GetIJ(0,0)*-1.0 +
            sides[s]->JTinv.GetIJ(0,1)*-1.0;
  }
  if (index==1)
  {

    value = sides[s]->JTinv.GetIJ(0,0)*1.0 +
            sides[s]->JTinv.GetIJ(0,1)*0.0;
  }

  value += beta*(sides[s]->JTinv.GetIJ(0,0)*0.0 +
                 sides[s]->JTinv.GetIJ(0,1)*1.0);


  return value;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GradVarphi_y
/**Precomputation of the partial derivative along y of the
 * shape function at a quadrature point.*/
double PolygonPWLFEValues::PreGradShape_y(int s, int i, int qpoint_index)
{
  int index = node_to_side_map[i][s];
  double value = 0;
  if (index==0)
  {

    value = sides[s]->JTinv.GetIJ(1,0)*-1.0 +
            sides[s]->JTinv.GetIJ(1,1)*-1.0;
  }
  if (index==1)
  {

    value = sides[s]->JTinv.GetIJ(1,0)*1.0 +
            sides[s]->JTinv.GetIJ(1,1)*0.0;
  }

  value += beta*(sides[s]->JTinv.GetIJ(1,0)*0.0 +
                 sides[s]->JTinv.GetIJ(1,1)*1.0);


  return value;
}


//###################################################################
/**Precomputes integrals of the shape functions.*/
void PolygonPWLFEValues::PreCompute()
{
  if (precomputed)
    return;

  // ==================================================== Precompute elements
  for (int s=0; s<num_of_subtris; s++)
  {
    for (int i=0; i<dofs; i++)
    {
      FEqp_data2d* pernode_data = new FEqp_data2d;
      for (int q=0; q<volume_quadrature.qpoints.size(); q++)
      {
        pernode_data->shape_qp.push_back(PreShape(s, i, q));
        pernode_data->gradshapex_qp.push_back(PreGradShape_x(s, i, q));
        pernode_data->gradshapey_qp.push_back(PreGradShape_y(s, i, q));
      }//for qp

      for (int q=0; q<surface_quadrature.qpoints.size(); q++)
      {
        //printf("%g\n",PreShape(s,i,q,ON_SURFACE));
        pernode_data->shape_qp_surf.push_back(PreShape(s,i,q,ON_SURFACE));
      }
      sides[s]->qp_data.push_back(pernode_data);
    }//for dof
  }//for side

  // ==================================================== Volume integrals
  for (int i=0; i<dofs; i++)
  {
    std::vector<double> gradijvalue_i(dofs, 0.0);
    std::vector<chi_mesh::Vector3> varphi_i_gradj(dofs, chi_mesh::Vector3());
    std::vector<double>           varphi_i_varphi_j(dofs,0);

    // Computing
    // GradVarphi_i*GradVarphi_j and
    // Varphi_i*GradVarphi_j_xyz and
    // Varphi_i*Varphi_j
    for (int j=0; j<dofs; j++)
    {
      gradijvalue_i[j] = 0.0;
      varphi_i_varphi_j[j] = 0.0;

      for (int s = 0; s < sides.size(); s++)
      {
        for (int qp=0; qp<volume_quadrature.qpoints.size();qp++)
        {
          gradijvalue_i[j]
            += volume_quadrature.weights[qp]*
               (GetGradShape_x(s, i, qp)*
                GetGradShape_x(s, j, qp) +
                GetGradShape_y(s, i, qp)*
                GetGradShape_y(s, j, qp))*
               DetJ(s,qp);
        }//for qp

        for (int qp=0; qp<volume_quadrature.qpoints.size();qp++)
        {
          double varphi_i = GetShape(s, i, qp);
          double weight = volume_quadrature.weights[qp];

          varphi_i_gradj[j].x
            += weight*varphi_i* GetGradShape_x(s, j, qp)*DetJ(s,qp);

          varphi_i_gradj[j].y
            += weight*varphi_i* GetGradShape_y(s, j, qp)*DetJ(s,qp);

          varphi_i_gradj[j].z = 0.0;

          varphi_i_varphi_j[j]
            += weight*varphi_i* GetShape(s, j, qp)*DetJ(s,qp);
        }// for qp
      }// for s
    }// for j

    //Computing Varphi_i
    double  valuei_i = 0.0;
    chi_mesh::Vector3 gradvalue_i(0, 0, 0);
    for (int s = 0; s < sides.size(); s++)
    {
      for (int qp=0; qp<volume_quadrature.qpoints.size();qp++)
      {
        valuei_i += volume_quadrature.weights[qp]*
                    GetShape(s, i, qp)*
                    DetJ(s,qp);
        gradvalue_i.x += volume_quadrature.weights[qp]*
                         GetGradShape_x(s, i, qp)*
                         DetJ(s,qp);
        gradvalue_i.y += volume_quadrature.weights[qp]*
                         GetGradShape_y(s, i, qp)*
                         DetJ(s,qp);
      }// for gp
    } // for s
    IntV_gradShapeI_gradShapeJ.push_back(gradijvalue_i);
    IntV_shapeI_gradshapeJ.push_back(varphi_i_gradj);
    IntV_shapeI_shapeJ.push_back(varphi_i_varphi_j);
    IntV_shapeI.push_back(valuei_i);
    IntV_gradshapeI.push_back(gradvalue_i);

    //=================================================== Surface integrals
    // Computing
    // Varphi_i*Varphi_j on each face and
    // Varphi_i on each face
    std::vector<std::vector<double>> varphi_i_varphi_j_surf;
    std::vector<std::vector<chi_mesh::Vector3>> varphi_i_gradvarphi_j_surf;
    std::vector<double> varphi_i_surf(num_of_subtris, 0.0);

    for (int f=0; f< num_of_subtris; f++)
    {
      std::vector<double>           f_varphi_i_varphi_j_surf(dofs,0);
      std::vector<chi_mesh::Vector3> f_varphi_i_grad_j_surf;
      f_varphi_i_grad_j_surf.resize(dofs,chi_mesh::Vector3());

      for (int j=0; j<dofs; j++)
      {
        double value_ij = 0.0;
        double value_x_ij = 0.0;
        double value_y_ij = 0.0;



        for (int qp=0; qp<surface_quadrature.qpoints.size();qp++)
        {
          value_ij
            += surface_quadrature.weights[qp]*
               GetShape(f, i, qp, ON_SURFACE)*
               GetShape(f, j, qp, ON_SURFACE)*
               DetJ(f,qp,ON_SURFACE);

          value_x_ij
            += surface_quadrature.weights[qp]*
               GetShape(f, i, qp, ON_SURFACE)*
               GetGradShape_x(f,j,qp)*
               DetJ(f,qp,ON_SURFACE);

          value_y_ij
            += surface_quadrature.weights[qp]*
               GetShape(f, i, qp, ON_SURFACE)*
               GetGradShape_y(f,j,qp)*
               DetJ(f,qp,ON_SURFACE);
        }// for gp

        f_varphi_i_varphi_j_surf[j] = value_ij;
        f_varphi_i_grad_j_surf[j].x = value_x_ij;
        f_varphi_i_grad_j_surf[j].y = value_y_ij;
        f_varphi_i_grad_j_surf[j].z = 0.0;
      }// for j
      varphi_i_varphi_j_surf.push_back(f_varphi_i_varphi_j_surf);
      varphi_i_gradvarphi_j_surf.push_back(f_varphi_i_grad_j_surf);

      double f_varphi_i_surf = 0.0;

      for (int qp=0; qp<surface_quadrature.qpoints.size();qp++)
      {
        f_varphi_i_surf
          += surface_quadrature.weights[qp]*
             GetShape(f, i, qp, ON_SURFACE)*
             DetJ(f,qp,ON_SURFACE);
      }// for gp

      varphi_i_surf[f] = f_varphi_i_surf;
    }// for f
    IntSi_shapeI_shapeJ.push_back(varphi_i_varphi_j_surf);
    IntSi_shapeI_gradshapeJ.push_back(varphi_i_gradvarphi_j_surf);
    IntS_shapeI.push_back(varphi_i_surf);
  }//for i



  //====================================== Reindexing surface integrals
//  IntS_shapeI_shapeJ.resize(num_of_subtris,
//                            std::vector<std::vector<double>>(dofs,std::vector<double>(dofs,0.0)));
  IntS_shapeI_shapeJ.resize(num_of_subtris);
  IntS_shapeI_gradshapeJ.resize(num_of_subtris);
  for (int f=0; f< num_of_subtris; f++)
  {
    IntS_shapeI_shapeJ[f].resize(dofs);
    IntS_shapeI_gradshapeJ[f].resize(dofs);
    for (int i=0; i<dofs; i++)
    {
      IntS_shapeI_shapeJ[f][i].resize(dofs);
      IntS_shapeI_gradshapeJ[f][i].resize(dofs);
      for (int j=0; j<dofs; j++)
      {
        IntS_shapeI_shapeJ[f][i][j] = IntSi_shapeI_shapeJ[i][f][j];
        IntS_shapeI_gradshapeJ[f][i][j] = IntSi_shapeI_gradshapeJ[i][f][j];
      }
    }
  }

//  double sum = 0.0;
//  double sum2= 0.0;
//  for (int i=0; i<dofs; i++)
//  {
//    for (int f=0; f<num_of_subtris; f++)
//    {
//      sum += IntS_shapeI[f][i];
//
//      for (int j=0; j<dofs; j++)
//        sum2 += IntS_shapeI_shapeJ[f][i][j];
//    }
//
//  }
//  printf("IntS %g\n",sum);
//  printf("IntS2 %g\n",sum2);
//  usleep(5000000);


  precomputed = true;
}