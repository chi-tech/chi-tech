#include "pwlc.h"

//###################################################################
/**Constructor.*/
SpatialDiscretization_PWLC::
  SpatialDiscretization_PWLC(chi_mesh::MeshContinuumPtr in_grid,
                             chi_math::finite_element::SetupFlags setup_flags/*=
                            chi_math::finite_element::SetupFlags::NO_FLAGS_SET*/,
                             chi_math::QuadratureOrder qorder/*=
                            chi_math::QuadratureOrder::SECOND*/) :
  SpatialDiscretization_FE(0, in_grid,
                           SDMType::PIECEWISE_LINEAR_CONTINUOUS,
                           setup_flags),
  line_quad_order_second(chi_math::QuadratureOrder::SECOND),
  tri_quad_order_second(chi_math::QuadratureOrder::SECOND),
  tet_quad_order_second(chi_math::QuadratureOrder::SECOND),
  line_quad_order_arbitrary(qorder),
  tri_quad_order_arbitrary (qorder),
  tet_quad_order_arbitrary (qorder)
{
  PreComputeCellSDValues();

  OrderNodes();
}
