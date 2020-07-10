#include "lbs_linear_boltzman_solver.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/** Computes the number of moments for the given mesher types*/
void LinearBoltzman::Solver::ComputeNumberOfMoments()
{
  if (options.geometry_type == GeometryType::ONED_SLAB)
  {
    int L = options.scattering_order;
    this->num_moments = L+1;
  }
  else if (options.geometry_type == GeometryType::TWOD_CARTESIAN or
           options.geometry_type == GeometryType::THREED_CARTESIAN)
  {
    int L = options.scattering_order;
    this->num_moments = L*(L+2) + 1;
  }
}

