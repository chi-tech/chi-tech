#include "spatial_discretization.h"

SpatialDiscretization::SpatialDiscretization(
  int in_dim,
  chi_mesh::MeshContinuumPtr in_grid,
  chi_math::SpatialDiscretizationType in_type) :
  dim(in_dim),
  type(in_type),
  ref_grid(in_grid)
{
}