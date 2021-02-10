#ifndef QUADRATURE_TETRAHEDRON_H
#define QUADRATURE_TETRAHEDRON_H

#include "quadrature.h"

namespace chi_math
{
  class QuadratureTetrahedron;
}

//###################################################################
/**Quadrature set for tetrahedrons.*/
class chi_math::QuadratureTetrahedron
{
public:
  std::vector<chi_math::QuadraturePointXYZ> qpoints;
  std::vector<double>     weights;

public:
  //00 Constructor
  explicit
  QuadratureTetrahedron(QuadratureOrder order,bool surface=false);

};

#endif