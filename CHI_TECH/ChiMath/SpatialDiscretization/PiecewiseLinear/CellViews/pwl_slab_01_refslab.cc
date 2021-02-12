#include "pwl_slab.h"

double SlabPWLFEView::SlabShape(int index, int qpoint_index, bool on_surface)
{
  double xi=0.0;
  if (!on_surface)
    xi = default_volume_quadrature.abscissae.at(qpoint_index);
  else
    if      (index == 0) xi=-1.0;
    else if (index == 1) xi= 1.0;

  double value = 0.0;
  if      (index == 0)
    value = 0.5*(1.0 - xi);
  else if (index == 1)
    value = 0.5*(1.0 + xi);

  return value;
}

double SlabPWLFEView::SlabGradShape(int index)
{
  double value = 0.0;

  if (index == 0)
    value = -1.0/h;
  else
    value =  1.0/h;

  return value;
}