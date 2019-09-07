#ifndef _quadrature_h
#define _quadrature_h

#include <stdio.h>
#include <vector>

namespace chi_math
{
  struct QPointXY;
  struct QPointXYZ;
  class CHI_QUADRATURE;
}

struct QPointXY
{
  double x;
  double y;
};

struct QPointXYZ
{
  double x;
  double y;
  double z;
};

//######################################################### Class def
/**Parent class for quadratures.*/
class CHI_QUADRATURE
{
public:
  std::vector<double> abscissae;
  std::vector<double> weights;

public:
  std::vector<double> q_weights;
  std::vector<double> q_points;
};


#endif