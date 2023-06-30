#include "pwl.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Get the number of ghost degrees-of-freedom.*/
size_t chi_math::SpatialDiscretization_PWLD::
GetNumGhostDOFs(const UnknownManager& unknown_manager) const
{
  return 0;
}

//###################################################################
/**Returns the ghost DOF indices.*/
std::vector<int64_t> chi_math::SpatialDiscretization_PWLD::
GetGhostDOFIndices(const UnknownManager& unknown_manager) const
{
  return {};
}

