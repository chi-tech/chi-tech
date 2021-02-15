#include "fv.h"

//###################################################################
/**Only constructor for this method.*/
SpatialDiscretization_FV::
  SpatialDiscretization_FV(chi_mesh::MeshContinuumPtr in_grid)
  : SpatialDiscretization(0, in_grid, SDMType::FINITE_VOLUME)
{
  mapping_initialized = false;
}