#include "pwl_polygon.h"

double PolygonMappingFE_PWL::TriShape(int index,
                                      const chi_mesh::Vector3& qpoint,
                                      bool on_surface/*false*/)
{
  double xi ;
  double eta;
  if (!on_surface)
  {
    xi = qpoint.x;
    eta= qpoint.y;
  }
  else
  {
    xi = 0.5*(qpoint[0] + 1.0);
    eta = 0.0;
  }

  double value = 0.0;
  if (index == 0)
    value = 1.0 - xi - eta;
  else if (index == 1)
    value = xi;
  else if (index == 2)
    value = eta;

  return value;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Varphi_x
/**Precomputation of the shape function at a quadrature point.*/
double PolygonMappingFE_PWL::SideShape(unsigned int side,
                                       unsigned int i,
                                       const chi_mesh::Vector3& qpoint,
                                       bool on_surface/*=false*/)
{
  int index = node_to_side_map[i][side];
  double value = 0.0;
  if (index == 0 or index == 1)
    value = TriShape(index,qpoint,on_surface);

  value += beta*TriShape(2,qpoint,on_surface);

  return value;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GradVarphi_x
/**Precomputation of the partial derivative along x of the
 * shape function at a quadrature point.*/
double PolygonMappingFE_PWL::SideGradShape_x(unsigned int side, int i)
{
  int index = node_to_side_map[i][side];
  double value = 0;
  if (index==0)
  {

    value = sides[side].JTinv.GetIJ(0, 0) * -1.0 +
            sides[side].JTinv.GetIJ(0, 1) * -1.0;
  }
  if (index==1)
  {

    value = sides[side].JTinv.GetIJ(0, 0) * 1.0 +
            sides[side].JTinv.GetIJ(0, 1) * 0.0;
  }

  value += beta*(sides[side].JTinv.GetIJ(0, 0) * 0.0 +
                 sides[side].JTinv.GetIJ(0, 1) * 1.0);


  return value;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GradVarphi_y
/**Precomputation of the partial derivative along y of the
 * shape function at a quadrature point.*/
double PolygonMappingFE_PWL::SideGradShape_y(unsigned int side, int i)
{
  int index = node_to_side_map[i][side];
  double value = 0;
  if (index==0)
  {

    value = sides[side].JTinv.GetIJ(1, 0) * -1.0 +
            sides[side].JTinv.GetIJ(1, 1) * -1.0;
  }
  if (index==1)
  {

    value = sides[side].JTinv.GetIJ(1, 0) * 1.0 +
            sides[side].JTinv.GetIJ(1, 1) * 0.0;
  }

  value += beta*(sides[side].JTinv.GetIJ(1, 0) * 0.0 +
                 sides[side].JTinv.GetIJ(1, 1) * 1.0);


  return value;
}