#include "quadrature_triangle.h"

#include <math.h>

//###################################################################
/**Initializes quadratures for use on triangles.*/
chi_math::QuadratureTriangle::QuadratureTriangle(int num_points, bool surface)
{
  if (surface)
  {
    if (num_points == 2)
    {
      abscissae.push_back( sqrt(3.0)/6.0 + 0.5);
      weights.push_back(0.5);

      abscissae.push_back(-sqrt(3.0)/6.0 + 0.5);
      weights.push_back(0.5);
    }
  }
  else
  {
    if (num_points == 3)
    {
      chi_math::QuadraturePointXY* new_qpoint;
      new_qpoint = new chi_math::QuadraturePointXY;
      new_qpoint->x = 1.0/6.0;
      new_qpoint->y = 1.0/6.0;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/6.0);

      new_qpoint = new chi_math::QuadraturePointXY;
      new_qpoint->x = 4.0/6.0;
      new_qpoint->y = 1.0/6.0;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/6.0);

      new_qpoint = new chi_math::QuadraturePointXY;
      new_qpoint->x = 1.0/6.0;
      new_qpoint->y = 4.0/6.0;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/6.0);
    }
  }
}