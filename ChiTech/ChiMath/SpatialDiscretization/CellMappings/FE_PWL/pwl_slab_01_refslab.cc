#include "pwl_slab.h"

double SlabMappingFE_PWL::SlabShape(int index, int qpoint_index, bool on_surface)
{
  double xi=0.0;
  if (!on_surface)
    xi = volume_quadrature.qpoints.at(qpoint_index)[0];
  else
    xi = static_cast<double>(index);

  double value = 0.0;
  if      (index == 0)
    value = 1.0 - xi;
  else if (index == 1)
    value = xi;

  return value;
}

double SlabMappingFE_PWL::SlabGradShape(int index)
{
  double value = 0.0;

  if (index == 0)
    value = -1.0/h;
  else
    value =  1.0/h;

  return value;
}