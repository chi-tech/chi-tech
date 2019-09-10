#include "quadrature_tetrahedron.h"

#include <math.h>

//###################################################################
/**Initialzes a set of points for a quadrature integration over
 * the volume of a tetrahedron.*/
chi_math::QuadratureTetrahedron::
 QuadratureTetrahedron(int num_points,bool surface)
{
  if (surface)
  {
    if (num_points == 3)
    {
      chi_math::QuadraturePointXYZ* new_qpoint;
      new_qpoint = new chi_math::QuadraturePointXYZ;
      new_qpoint->x = 1.0/6.0;
      new_qpoint->y = 1.0/6.0;
      new_qpoint->z = 0.0;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/6.0);

      new_qpoint = new chi_math::QuadraturePointXYZ;
      new_qpoint->x = 4.0/6.0;
      new_qpoint->y = 1.0/6.0;
      new_qpoint->z = 0.0;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/6.0);

      new_qpoint = new chi_math::QuadraturePointXYZ;
      new_qpoint->x = 1.0/6.0;
      new_qpoint->y = 4.0/6.0;
      new_qpoint->z = 0.0;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/6.0);
    }
  }
  else
  {
    if (num_points == 1)
    {
      chi_math::QuadraturePointXYZ* new_qpoint;
      new_qpoint = new chi_math::QuadraturePointXYZ;
      new_qpoint->x = 0.25;
      new_qpoint->y = 0.25;
      new_qpoint->z = 0.25;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/6);
    }
    else if (num_points == 3)
    {
      double a = (5.0+3.0*sqrt(5))/20.0;
      double b = (5.0-    sqrt(5))/20.0;
      chi_math::QuadraturePointXYZ* new_qpoint;
      new_qpoint = new chi_math::QuadraturePointXYZ;
      new_qpoint->x = a;
      new_qpoint->y = b;
      new_qpoint->z = b;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/24.0);

      new_qpoint = new chi_math::QuadraturePointXYZ;
      new_qpoint->x = b;
      new_qpoint->y = a;
      new_qpoint->z = b;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/24.0);

      new_qpoint = new chi_math::QuadraturePointXYZ;
      new_qpoint->x = b;
      new_qpoint->y = b;
      new_qpoint->z = a;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/24.0);

      new_qpoint = new chi_math::QuadraturePointXYZ;
      new_qpoint->x = b;
      new_qpoint->y = b;
      new_qpoint->z = b;
      qpoints.push_back(new_qpoint);
      weights.push_back(1.0/24.0);
    }
  }

}