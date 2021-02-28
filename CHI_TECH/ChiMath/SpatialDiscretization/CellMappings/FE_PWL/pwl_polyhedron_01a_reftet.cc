#include "pwl_polyhedron.h"

double PolyhedronMappingFE_PWL::TetShape(unsigned int index,
                                         const chi_mesh::Vector3& qpoint,
                                         bool on_surface/*=false*/)
{
  double value = 0.0;

  if (index == 0)
  {value = 1.0 - qpoint.x - qpoint.y - qpoint.z;}
  if (index == 1){value = qpoint.x;}
  if (index == 2){value = qpoint.y;}
  if (index == 3){value = qpoint.z;}

  return value;
}

double PolyhedronMappingFE_PWL::TetGradShape_x(const unsigned int index)
{
  double value = 0.0;
  if (index == 0){value = -1.0;}
  if (index == 1){value =  1.0;}
  if (index == 2){value =  0.0;}
  if (index == 3){value =  0.0;}

  return value;
}

double PolyhedronMappingFE_PWL::TetGradShape_y(const unsigned int index)
{
  double value = 0.0;
  if (index == 0){value = -1.0;}
  if (index == 1){value =  0.0;}
  if (index == 2){value =  1.0;}
  if (index == 3){value =  0.0;}

  return value;
}

double PolyhedronMappingFE_PWL::TetGradShape_z(const unsigned int index)
{
  double value = 0.0;
  if (index == 0){value = -1.0;}
  if (index == 1){value =  0.0;}
  if (index == 2){value =  0.0;}
  if (index == 3){value =  1.0;}

  return value;
}