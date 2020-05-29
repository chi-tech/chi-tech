#include "fv.h"

//###################################################################
/**Only constructor for this method.*/
SpatialDiscretization_FV::SpatialDiscretization_FV(int dim)
  : SpatialDiscretization(dim)
{
  mapping_initialized = false;
}