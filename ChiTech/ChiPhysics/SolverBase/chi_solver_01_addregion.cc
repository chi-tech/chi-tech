#include "chi_solver.h"

//#########################################################
/**Adds a region to a solver.*/
void chi_physics::Solver::AddRegion(chi_mesh::Region* region)
{
  this->regions.push_back(region);
}
