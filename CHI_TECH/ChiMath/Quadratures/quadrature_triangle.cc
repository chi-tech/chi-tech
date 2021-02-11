#include "quadrature_triangle.h"

#include <cmath>
#include <stdexcept>

//###################################################################
/**Initializes quadratures for use on triangles.*/
chi_math::QuadratureTriangle::QuadratureTriangle(QuadratureOrder order, bool surface)
{
  double x=0.0,y=0.0,z=0.0;
  if (surface)
  {
    switch (order)
    {
      case QuadratureOrder::CONSTANT:
      case QuadratureOrder::FIRST:
      {
        qpoints.emplace_back(0.5,0.0,0.0);
        weights.emplace_back(1.0);
        break;
      }
      case QuadratureOrder::SECOND:
      {
        qpoints.emplace_back( sqrt(3.0)/6.0 + 0.5,0.0,0.0);
        weights.push_back(0.5);

        qpoints.emplace_back(-sqrt(3.0)/6.0 + 0.5,0.0,0.0);
        weights.push_back(0.5);
        break;
      }
      default:
        throw std::invalid_argument("Invalid order specified for surface "
                                    "QuadratureTriangle.");
    }//switch order
  }
  else
  {
    switch (order)
    {
      case QuadratureOrder::CONSTANT:
      case QuadratureOrder::FIRST:
      {
        x = 1.0/3.0;
        y = 1.0/3.0;
        qpoints.emplace_back(x,y,z);
        weights.push_back(0.5);
        break;
      }
      case QuadratureOrder::SECOND:
      {
        x = 1.0/6.0;
        y = 1.0/6.0;
        qpoints.emplace_back(x,y,z);
        weights.push_back(1.0/6.0);

        x = 4.0/6.0;
        y = 1.0/6.0;
        qpoints.emplace_back(x,y,z);
        weights.push_back(1.0/6.0);

        x = 1.0/6.0;
        y = 4.0/6.0;
        qpoints.emplace_back(x,y,z);
        weights.push_back(1.0/6.0);
        break;
      }
      default:
        throw std::invalid_argument("Invalid order specified for "
                                    "QuadratureTriangle.");
    }//switch order
  }
}