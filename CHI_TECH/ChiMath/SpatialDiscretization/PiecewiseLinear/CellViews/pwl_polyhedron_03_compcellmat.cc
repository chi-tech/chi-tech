#include "pwl_polyhedron.h"

/**Precomputes the shape function values of a face-side pair
 * at a quadrature point*/
double PolyhedronPWLFEValues::PreShape(int face_index, int side_index,
                                       int i, int qpoint_index, bool on_surface)
{
  double value = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  value += TetShape(index, qpoint_index, on_surface);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    value += betaf* TetShape(1, qpoint_index, on_surface);}
  value += alphac* TetShape(3, qpoint_index, on_surface);

  return value;
}

/**Precomputes the gradx-shape function values of a face-side pair
 * at a quadrature point*/
double PolyhedronPWLFEValues::PreGradShape_x(int face_index,
                                             int side_index,
                                             int i)
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdx += betaf* TetGradShape_x(1);}
  tetdfdx += alphac* TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdy += betaf* TetGradShape_y(1);}
  tetdfdy += alphac* TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdz += betaf* TetGradShape_z(1);}
  tetdfdz += alphac* TetGradShape_z(3);

  value += face_data[face_index].sides[side_index].JTinv.GetIJ(0, 0) * tetdfdx;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(0, 1) * tetdfdy;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(0, 2) * tetdfdz;

  return value;
}

/**Precomputes the grady-shape function values of a face-side pair
 * at a quadrature point*/
double PolyhedronPWLFEValues::PreGradShape_y(int face_index,
                                             int side_index,
                                             int i)
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdx += betaf* TetGradShape_x(1);}
  tetdfdx += alphac* TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdy += betaf* TetGradShape_y(1);}
  tetdfdy += alphac* TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdz += betaf* TetGradShape_z(1);}
  tetdfdz += alphac* TetGradShape_z(3);

  value += face_data[face_index].sides[side_index].JTinv.GetIJ(1, 0) * tetdfdx;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(1, 1) * tetdfdy;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(1, 2) * tetdfdz;

  return value;
}

/**Precomputes the gradz-shape function values of a face-side pair
 * at a quadrature point*/
double PolyhedronPWLFEValues::PreGradShape_z(int face_index,
                                             int side_index,
                                             int i)
{
  double value = 0.0;
  double tetdfdx = 0.0;
  double tetdfdy = 0.0;
  double tetdfdz = 0.0;
  int    index = node_side_maps[i].face_map[face_index].
    side_map[side_index].index;
  double betaf = face_betaf[face_index];

  tetdfdx += TetGradShape_x(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdx += betaf* TetGradShape_x(1);}
  tetdfdx += alphac* TetGradShape_x(3);

  tetdfdy += TetGradShape_y(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdy += betaf* TetGradShape_y(1);}
  tetdfdy += alphac* TetGradShape_y(3);

  tetdfdz += TetGradShape_z(index);
  if (node_side_maps[i].face_map[face_index].side_map[side_index].part_of_face){
    tetdfdz += betaf* TetGradShape_z(1);}
  tetdfdz += alphac* TetGradShape_z(3);

  value += face_data[face_index].sides[side_index].JTinv.GetIJ(2, 0) * tetdfdx;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(2, 1) * tetdfdy;
  value += face_data[face_index].sides[side_index].JTinv.GetIJ(2, 2) * tetdfdz;

  return value;
}




/**Precomputes cell volume and surface integrals.*/
void PolyhedronPWLFEValues::PreCompute()
{
  // ==================================================== Precompute elements
  // Precomputing the values of shape functions and
  // derivatives of shape functions at quadrature points
  // for each tetrahedron
  if (precomputed){return; }

  std::vector<std::vector<std::vector<double>>> IntSi_shapeI_shapeJ;
  std::vector<std::vector<std::vector<chi_mesh::Vector3>>> IntSi_shapeI_gradshapeJ;

  for (size_t f=0; f < face_data.size(); f++)
  {
    for (size_t s=0; s < face_data[f].sides.size(); s++)
    {
      for (size_t i=0; i<dofs; i++)
      {
        FEqp_data3d pernode_data;

        //========== Reserving data
        pernode_data.gradshapex_qp.reserve(
          quadratures[DEG3]->qpoints.size());
        pernode_data.gradshapey_qp.reserve(
          quadratures[DEG3]->qpoints.size());
        pernode_data.gradshapez_qp.reserve(
          quadratures[DEG3]->qpoints.size());
        pernode_data.shape_qp.reserve(
          quadratures[DEG3]->qpoints.size());
        pernode_data.shape_qp_surf.reserve(
          quadratures[DEG3_SURFACE]->qpoints.size());

        //Prestore GradVarphi_xyz
        for (size_t qp=0; qp<quadratures[DEG3]->qpoints.size(); qp++)
        {
          pernode_data.gradshapex_qp.push_back(PreGradShape_x(f, s, i));
          pernode_data.gradshapey_qp.push_back(PreGradShape_y(f, s, i));
          pernode_data.gradshapez_qp.push_back(PreGradShape_z(f, s, i));
        }
        //Prestore Varphi
        for (size_t qp=0; qp<quadratures[DEG3]->qpoints.size(); qp++)
        {
          pernode_data.shape_qp.push_back(PreShape(f, s, i, qp));
        }
        //Prestore Varphi on surface
        for (size_t qp=0; qp<quadratures[DEG3_SURFACE]->qpoints.size(); qp++)
        {
          pernode_data.shape_qp_surf.push_back(PreShape(f, s, i, qp, ON_SURFACE));
        }


        face_data[f].sides[s].qp_data.push_back(pernode_data);
      } // for i

    } //for side
  } //for face

  // ==================================================== Volume integrals
  IntV_gradShapeI_gradShapeJ.reserve(dofs);
  IntV_shapeI_gradshapeJ.reserve(dofs);
  IntV_shapeI_shapeJ.reserve(dofs);
  IntV_shapeI.reserve(dofs);

//  IntV_gradShapeI_gradShapeJ.push_back(gradijvalue_i);
//  IntV_shapeI_gradshapeJ.push_back(varphi_i_gradj);
//  IntV_shapeI_shapeJ.push_back(varphi_i_varphi_j);
//  IntV_shapeI.push_back(valuei_i);
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

      for (size_t f=0; f < face_data.size(); f++)
      {
        for (size_t s = 0; s < face_data[f].sides.size(); s++)
        {
          for (size_t qp=0; qp<quadratures[DEG3]->qpoints.size();qp++)
          {
            gradijvalue_i[j]
              += quadratures[DEG3]->weights[qp]*
                 (GetGradShape_x(f, s, i, qp)*GetGradShape_x(f, s, j, qp) +
                  GetGradShape_y(f, s, i, qp)*GetGradShape_y(f, s, j, qp) +
                  GetGradShape_z(f, s, i, qp)*GetGradShape_z(f, s, j, qp))*
                 DetJ(f,s,qp);
          }//for qp

          for (size_t qp=0; qp<quadratures[DEG3]->qpoints.size();qp++)
          {
            double varphi_i = GetShape(f, s, i, qp);
            double weight = quadratures[DEG3]->weights[qp];

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
        for (size_t qp=0; qp<quadratures[DEG3]->qpoints.size();qp++)
        {
          valuei_i += quadratures[DEG3]->weights[qp]*
            GetShape(f, s, i, qp)*
                      DetJ(f,s,qp);
        }// for gp
      } // for s
    }// for f
    IntV_gradShapeI_gradShapeJ.push_back(std::move(gradijvalue_i));
    IntV_shapeI_gradshapeJ.push_back(std::move(varphi_i_gradj));
    IntV_shapeI_shapeJ.push_back(std::move(varphi_i_varphi_j));
    IntV_shapeI.push_back(std::move(valuei_i));

    //=================================================== Surface integrals
    // Computing
    // Varphi_i*Varphi_j on each face and
    // Varphi_i on each face
    std::vector<std::vector<double>> varphi_i_varphi_j_surf;
    std::vector<std::vector<chi_mesh::Vector3>> varphi_i_gradvarphi_j_surf;
    std::vector<double> varphi_i_surf(face_data.size(), 0.0);

    for (size_t f=0; f < face_data.size(); f++)
    {
      std::vector<double> f_varphi_i_varphi_j_surf(dofs,0);
      std::vector<chi_mesh::Vector3> f_varphi_i_grad_j_surf;
      f_varphi_i_grad_j_surf.resize(dofs,chi_mesh::Vector3());

      for (int j=0; j<dofs; j++)
      {
        double value_ij = 0.0;
        double value_x_ij = 0.0;
        double value_y_ij = 0.0;
        double value_z_ij = 0.0;

        for (size_t s = 0; s < face_data[f].sides.size(); s++)
        {
          for (size_t qp=0; qp<quadratures[DEG3_SURFACE]->qpoints.size();qp++)
          {
            value_ij
              += quadratures[DEG3_SURFACE]->weights[qp]*
              GetShape(f, s, i, qp, ON_SURFACE)*
              GetShape(f, s, j, qp, ON_SURFACE)*
                 DetJ(f,s,qp,ON_SURFACE);

            value_x_ij
              += quadratures[DEG3_SURFACE]->weights[qp]*
                 GetShape(f, s, i, qp, ON_SURFACE)*
                 GetGradShape_x(f, s, j, qp)*
                 DetJ(f,s,qp,ON_SURFACE);

            value_y_ij
              += quadratures[DEG3_SURFACE]->weights[qp]*
                 GetShape(f, s, i, qp, ON_SURFACE)*
                 GetGradShape_y(f, s, j, qp)*
                 DetJ(f,s,qp,ON_SURFACE);

            value_z_ij
              += quadratures[DEG3_SURFACE]->weights[qp]*
                 GetShape(f, s, i, qp, ON_SURFACE)*
                 GetGradShape_z(f, s, j, qp)*
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
        for (size_t qp=0; qp<quadratures[DEG3_SURFACE]->qpoints.size();qp++)
        {
          f_varphi_i_surf
            += quadratures[DEG3_SURFACE]->weights[qp]*
            GetShape(f, s, i, qp, ON_SURFACE)*
               DetJ(f,s,qp,ON_SURFACE);
        }// for gp
      }
      varphi_i_surf[f] = f_varphi_i_surf;
    }// for f
    IntSi_shapeI_shapeJ.push_back(std::move(varphi_i_varphi_j_surf));
    IntSi_shapeI_gradshapeJ.push_back(std::move(varphi_i_gradvarphi_j_surf));
    IntS_shapeI.push_back(std::move(varphi_i_surf));
  }// for i

  //====================================== Reindexing surface integrals
  IntS_shapeI_shapeJ.resize(face_data.size());
  IntS_shapeI_gradshapeJ.resize(face_data.size());
  for (size_t f=0; f < face_data.size(); f++)
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

  precomputed = true;
}