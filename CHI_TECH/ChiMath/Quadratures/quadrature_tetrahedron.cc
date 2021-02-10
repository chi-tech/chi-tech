#include "quadrature_tetrahedron.h"

#include <cmath>
#include <stdexcept>

//###################################################################
/**Initialzes a set of points for a quadrature integration over
 * the volume of a tetrahedron.*/
chi_math::QuadratureTetrahedron::
 QuadratureTetrahedron(QuadratureOrder order,bool surface)
{
  double x=0.0,y=0.0,z=0.0;

  if (surface)
  {
    switch (order)
    {
      case QuadratureOrder::FIRST:
      {
        x = 1.0/3.0;
        y = 1.0/3.0;
        z = 0.0;
        qpoints.emplace_back(x,y,z);
        weights.push_back(0.5);
        break;
      }
      case QuadratureOrder::SECOND:
      {
        x = 1.0/6.0;
        y = 1.0/6.0;
        z = 0.0;
        qpoints.emplace_back(x,y,z);
        weights.push_back(1.0/6.0);

        x = 4.0/6.0;
        y = 1.0/6.0;
        z = 0.0;
        qpoints.emplace_back(x,y,z);
        weights.push_back(1.0/6.0);

        x = 1.0/6.0;
        y = 4.0/6.0;
        z = 0.0;
        qpoints.emplace_back(x,y,z);
        weights.push_back(1.0/6.0);
        break;
      }
      default:
      {
        throw std::invalid_argument("Invalid surface order specified for "
                                    "QuadratureTetrahedron.");
      }
    }//switch order
  }
  else
  {
    switch (order)
    {
      case QuadratureOrder::FIRST:
      {
        x = 0.25;
        y = 0.25;
        z = 0.25;
        qpoints.emplace_back(x,y,z);
        weights.push_back(1.0/6.0);
        break;
      }
      case QuadratureOrder::SECOND:
      {
        double a = (5.0+3.0*sqrt(5))/20.0;
        double b = (5.0-    sqrt(5))/20.0;

        x = a;
        y = b;
        z = b;
        qpoints.emplace_back(x,y,z);
        weights.push_back(1.0/24.0);

        x = b;
        y = a;
        z = b;
        qpoints.emplace_back(x,y,z);
        weights.push_back(1.0/24.0);

        x = b;
        y = b;
        z = a;
        qpoints.emplace_back(x,y,z);
        weights.push_back(1.0/24.0);

        x = b;
        y = b;
        z = b;
        qpoints.emplace_back(x,y,z);
        weights.push_back(1.0/24.0);
        break;
      }
      default:
      {
        throw std::invalid_argument("Invalid order specified for "
                                    "QuadratureTetrahedron.");
      }
    }//switch order
  }//not surface

}