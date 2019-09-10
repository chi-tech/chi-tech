#ifndef _chi_diffusion_bndry_h
#define _chi_diffusion_bndry_h

#define PROCESS_BOUNDARY 1
#define NO_BOUNDARY           0
#define DIFFUSION_REFLECTING -1
#define DIFFUSION_DIRICHLET  -2
#define DIFFUSION_NEUMANN    -3
#define DIFFUSION_VACUUM     -4
#define DIFFUSION_ROBIN      -5

#include "../chi_diffusion.h"

//###################################################################
/**Parent class for diffusion boundaries*/
class chi_diffusion::Boundary
{
public:
  int type;

};


#endif