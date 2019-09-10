#ifndef _chi_diffusion_bndry_reflecting_h
#define _chi_diffusion_bndry_reflecting_h

#include "chi_diffusion_bndry.h"

//###################################################################
/**Reflecting boundary condition.*/
class chi_diffusion::BoundaryReflecting : public chi_diffusion::Boundary
{
public:
  BoundaryReflecting()
  {
    type = DIFFUSION_REFLECTING;
  }
};


#endif