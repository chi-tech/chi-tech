#include "pwl_base.h"

//###################################################################
/**Constructor*/
chi_math::SpatialDiscretization_PWLBase::
  SpatialDiscretization_PWLBase(const chi_mesh::MeshContinuum& in_grid,
                                finite_element::SetupFlags in_setup_flags,
                                QuadratureOrder qorder,
                                SDMType in_type,
                                CoordinateSystemType in_cs_type
                                ) :
    SpatialDiscretization_FE(in_grid, in_cs_type, in_type,in_setup_flags),
    line_quad_order_arbitrary_(qorder),
    tri_quad_order_arbitrary_(qorder),
    quad_quad_order_arbitrary_(qorder),
    tet_quad_order_arbitrary_(qorder)
{}