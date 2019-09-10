#ifndef _quadrature_triangle_h
#define _quadrature_triangle_h

#include "quadrature.h"

namespace chi_math
{
  class QuadratureTriangle;
}

class chi_math::QuadratureTriangle : public chi_math::Quadrature
{
public:
  std::vector<chi_math::QuadraturePointXY*>  qpoints;
  std::vector<double>     weights;

public:
  QuadratureTriangle(int num_points=3, bool surface=false);
};

#endif