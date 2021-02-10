#ifndef QUADRATURE_TRIANGLE_H
#define QUADRATURE_TRIANGLE_H

#include "quadrature.h"

namespace chi_math
{
  class QuadratureTriangle;
}

class chi_math::QuadratureTriangle
{
public:
  std::vector<chi_math::QuadraturePointXYZ>  qpoints;
  std::vector<double>     weights;

public:
  explicit
  QuadratureTriangle(QuadratureOrder order, bool surface=false);
};

#endif