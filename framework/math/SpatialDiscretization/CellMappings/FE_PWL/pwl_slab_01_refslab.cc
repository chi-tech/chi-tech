#include "pwl_slab.h"

double chi_math::SlabMappingFE_PWL::SlabShape(int index,
                                              const chi_mesh::Vector3& qpoint,
                                              bool on_surface,
                                              const int edge) const
{
  double xi=0.0;
  if (!on_surface)
    xi = qpoint.x;
  else
    xi = static_cast<double>(edge);

  double value = 0.0;
  if      (index == 0)
    value = 1.0 - xi;
  else if (index == 1)
    value = xi;

  return value;
}

double chi_math::SlabMappingFE_PWL::SlabGradShape(int index) const
{
  double value = 0.0;

  if (index == 0)
    value = -1.0 / h_;
  else
    value = 1.0 / h_;

  return value;
}