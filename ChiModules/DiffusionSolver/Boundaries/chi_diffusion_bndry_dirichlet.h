#ifndef _chi_diffusion_bndry_dirichlet_h
#define _chi_diffusion_bndry_dirichlet_h

#include "chi_diffusion_bndry.h"

//###################################################################
/**Dirichlet boundary.*/
class chi_diffusion::BoundaryDirichlet : public chi_diffusion::Boundary
{
public:
  double boundary_value;

public:
  BoundaryDirichlet() {
    type = DIFFUSION_DIRICHLET;
    boundary_value = 0.0;
  }
};

#endif