#include "pwl.h"

//###################################################################
/**Constructor.*/
SpatialDiscretization_PWL::
  SpatialDiscretization_PWL(int dim,
                            chi_math::SpatialDiscretizationType sd_method) :
  SpatialDiscretization(dim, sd_method),
  line_quad_order_second(chi_math::QuadratureOrder::SECOND),
  tri_quad_order_second(chi_math::QuadratureOrder::SECOND),
  tet_quad_order_second(chi_math::QuadratureOrder::SECOND)
{
  mapping_initialized = false;
}

