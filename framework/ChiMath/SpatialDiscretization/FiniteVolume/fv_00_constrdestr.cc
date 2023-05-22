#include "fv.h"

//###################################################################
/**Only constructor for this method.*/
chi_math::SpatialDiscretization_FV::
  SpatialDiscretization_FV(const chi_mesh::MeshContinuum& in_grid,
                           chi_math::CoordinateSystemType in_cs_type) :
  chi_math::SpatialDiscretization(in_grid, in_cs_type, SDMType::FINITE_VOLUME)
{
  CreateCellMappings();

  OrderNodes();
}