#ifndef MG_DIFFUSION_BOUNDARY_H
#define MG_DIFFUSION_BOUNDARY_H

#include "array"

namespace mg_diffusion
{
  class Boundary;
  
  enum class BoundaryType : int
  {
    Reflecting = 1,
    Neumann    = 3,
    Robin      = 4,
    Vacuum     = 5
  };
}

//###################################################################
/**Parent class for multigroup diffusion boundaries*/
class mg_diffusion::Boundary
{
  public :
  BoundaryType type = BoundaryType::Vacuum;

  // std::array<std::vector<double>, 3> mg_values2;
  std::array<double, 3> mg_values = {0.,0.,0.};

};

#endif //MG_DIFFUSION_BOUNDARY_H