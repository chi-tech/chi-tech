#include "pwl_polygon.h"

double PolygonPWLFEValues::TriShape(int index, int qpoint_index,
                                    bool on_surface/*false*/)
{
  double xi  = 0.0;
  double eta = 0.0;
  if (!on_surface)
  {
    auto& qpoint = default_volume_quadrature.qpoints.at(qpoint_index);

    xi = qpoint.x;
    eta= qpoint.y;
  }
  else
  {
    xi = 0.5*(default_surface_quadrature.qpoints[qpoint_index][0] + 1.0);
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
double PolygonPWLFEValues::SideShape(int side, int i, int qpoint_index,
                                     bool on_surface/*=false*/)
{
  int index = node_to_side_map[i][side];
  double value = 0.0;
  if (index == 0 or index == 1)
    value = TriShape(index,qpoint_index,on_surface);

  value += beta*TriShape(2,qpoint_index,on_surface);

  return value;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GradVarphi_x
/**Precomputation of the partial derivative along x of the
 * shape function at a quadrature point.*/
double PolygonPWLFEValues::SideGradShape_x(int side, int i)
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
double PolygonPWLFEValues::SideGradShape_y(int side, int i)
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