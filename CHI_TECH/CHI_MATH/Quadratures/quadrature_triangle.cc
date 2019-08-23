#include "quadrature_triangle.h"

#include <math.h>

//###################################################################
/**Initializes quadratures for use on triangles.*/
CHI_QUADRATURE_TRIANGLE::CHI_QUADRATURE_TRIANGLE(int num_points, bool surface)
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
      QPointXY* new_qpoint;
      new_qpoint = new QPointXY;
      new_qpoint->x = 1.0/6.0;
      new_qpoint->y = 1.0/6.0;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/6.0);

      new_qpoint = new QPointXY;
      new_qpoint->x = 4.0/6.0;
      new_qpoint->y = 1.0/6.0;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/6.0);

      new_qpoint = new QPointXY;
      new_qpoint->x = 1.0/6.0;
      new_qpoint->y = 4.0/6.0;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/6.0);
    }
  }
}