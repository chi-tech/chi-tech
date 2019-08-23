#ifndef _quadrature_tetrahedron_h
#define _quadrature_tetrahedron_h

#include "quadrature.h"

//###################################################################
/**Quadrature set for tetrahedrons.*/
class CHI_QUADRATURE_TETRAHEDRON : public CHI_QUADRATURE
{
public:
  std::vector<QPointXYZ*> qpoints;
  std::vector<double>     weights;

public:
  //00 Constructor
  CHI_QUADRATURE_TETRAHEDRON(int num_points=4,bool surface=false);

};

#endif