#include "spatial_discretization.h"

SpatialDiscretization::SpatialDiscretization(
  int dim,
  chi_math::SpatialDiscretizationType in_type) :
 type(in_type)
{
  this->dim = dim;
}