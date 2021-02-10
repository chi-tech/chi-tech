#include "pwl.h"

//###################################################################
/**Constructor.*/
SpatialDiscretization_PWL::
  SpatialDiscretization_PWL(int dim,
                            chi_math::SpatialDiscretizationType sd_method) :
  SpatialDiscretization(dim, sd_method),
  tri_quad_deg5(chi_math::QuadratureOrder::SECOND),
  tri_quad_deg3_surf(chi_math::QuadratureOrder::SECOND,true),
  tet_quad_order2(chi_math::QuadratureOrder::SECOND),
  tet_quad_order2_surface(chi_math::QuadratureOrder::SECOND,true)
{
  mapping_initialized = false;
}

