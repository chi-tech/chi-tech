#ifndef _quadrature_tetrahedron_h
#define _quadrature_tetrahedron_h

#include "quadrature.h"

namespace chi_math
{
  class QuadratureTetrahedron;
}

//###################################################################
/**Quadrature set for tetrahedrons.*/
class chi_math::QuadratureTetrahedron : public chi_math::Quadrature
{
public:
  std::vector<chi_math::QuadraturePointXYZ*> qpoints;
  std::vector<double>     weights;

public:
  //00 Constructor
  QuadratureTetrahedron(int num_points=4,bool surface=false);

};

#endif