#include "fv.h"

//###################################################################
/**Only constructor for this method.*/
chi_math::SpatialDiscretization_FV::
  SpatialDiscretization_FV(chi_mesh::MeshContinuumPtr& in_grid,
                           chi_math::CoordinateSystemType in_cs_type) :
  chi_math::SpatialDiscretization(in_grid, in_cs_type, SDMType::FINITE_VOLUME)
{
  CreateCellMappings();

  OrderNodes();
}