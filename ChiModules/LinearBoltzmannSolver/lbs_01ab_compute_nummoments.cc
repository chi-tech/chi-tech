#include "lbs_linear_boltzmann_solver.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/** Computes the number of moments for the given mesher types*/
void LinearBoltzmann::Solver::ComputeNumberOfMoments()
{
  int L = options.scattering_order;

  if (options.geometry_type == GeometryType::ONED_SLAB)
  {
    this->num_moments = L+1;
  }
  else if (options.geometry_type == GeometryType::TWOD_CARTESIAN)
  {
    this->num_moments = 0;
    for (int ell=0; ell<=L; ell++)
      for (int m=-ell; m<=ell; m+=2)
        if (ell == 0 or m != 0)
          ++this->num_moments;
  }
  else if (options.geometry_type == GeometryType::THREED_CARTESIAN)
  {
    this->num_moments = (L+1)*(L+1);
  }
}

