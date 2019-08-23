#ifndef _quadrature_h
#define _quadrature_h

#include <stdio.h>
#include <vector>

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