#ifndef _quadrature_triangle_h
#define _quadrature_triangle_h

#include "quadrature.h"

class CHI_QUADRATURE_TRIANGLE : public CHI_QUADRATURE
{
public:
  std::vector<QPointXY*>  qpoints;
  std::vector<double>     weights;

public:
  CHI_QUADRATURE_TRIANGLE(int num_points=3, bool surface=false);
};

#endif