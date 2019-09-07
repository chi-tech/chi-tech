#ifndef _quadrature_h
#define _quadrature_h

#include <stdio.h>
#include <vector>

namespace chi_math
{
  struct QuadraturePointXY;
  struct QuadraturePointXYZ;
  class Quadrature;
}

struct chi_math::QuadraturePointXY
{
  double x;
  double y;
};

struct chi_math::QuadraturePointXYZ
{
  double x;
  double y;
  double z;
};

//######################################################### Class def
/**Parent class for quadratures.*/
class chi_math::Quadrature
{
public:
  std::vector<double> abscissae;
  std::vector<double> weights;

public:
  std::vector<double> q_weights;
  std::vector<double> q_points;
};


#endif