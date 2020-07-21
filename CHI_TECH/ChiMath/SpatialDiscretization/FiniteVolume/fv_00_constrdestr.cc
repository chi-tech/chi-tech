#include "fv.h"

//###################################################################
/**Only constructor for this method.*/
SpatialDiscretization_FV::
  SpatialDiscretization_FV(int dim,
                           chi_math::SpatialDiscretizationType sd_method)
  : SpatialDiscretization(dim, sd_method)
{
  mapping_initialized = false;
}