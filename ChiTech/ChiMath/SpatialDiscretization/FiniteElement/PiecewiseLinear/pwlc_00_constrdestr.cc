#include "pwlc.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Constructor.*/
SpatialDiscretization_PWLC::
  SpatialDiscretization_PWLC(chi_mesh::MeshContinuumPtr& in_grid,
                             chi_math::finite_element::SetupFlags setup_flags,
                             chi_math::QuadratureOrder qorder,
                             chi_math::CoordinateSystemType in_cs_type) :
  SpatialDiscretization_FE(0, in_grid, in_cs_type,
                           SDMType::PIECEWISE_LINEAR_CONTINUOUS,
                           setup_flags),
  line_quad_order_arbitrary(qorder),
  tri_quad_order_arbitrary(qorder),
  quad_quad_order_arbitrary(qorder),
  tet_quad_order_arbitrary(qorder),
  hex_quad_order_arbitrary(qorder)
{
  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Creating Piecewise Linear Continuous "
                   "Finite Element spatial discretizaiton.";

  if (setup_flags == chi_math::finite_element::COMPUTE_UNIT_INTEGRALS)
  {
    const auto qorder_min = (int)chi_math::QuadratureOrder::SECOND;

    if ((int)line_quad_order_arbitrary.order < qorder_min)
      chi_log.Log(LOG_ALLWARNING)
        << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
        << "(int)line_quad_order_arbitrary.order < " << qorder_min << ".";

    if ((int)tri_quad_order_arbitrary.order < qorder_min)
      chi_log.Log(LOG_ALLWARNING)
        << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
        << "(int)tri_quad_order_arbitrary.order < " << qorder_min << ".";

    if ((int)quad_quad_order_arbitrary.order < qorder_min)
      chi_log.Log(LOG_ALLWARNING)
        << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
        << "(int)quad_quad_order_arbitrary.order < " << qorder_min << ".";

    if ((int)tet_quad_order_arbitrary.order < qorder_min)
      chi_log.Log(LOG_ALLWARNING)
        << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
        << "(int)tet_quad_order_arbitrary.order < " << qorder_min << ".";

    if ((int)hex_quad_order_arbitrary.order < qorder_min)
      chi_log.Log(LOG_ALLWARNING)
        << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
        << "(int)hex_quad_order_arbitrary.order < " << qorder_min << ".";
  }

  if (setup_flags != chi_math::finite_element::NO_FLAGS_SET)
    PreComputeCellSDValues();

  OrderNodes();
  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Done creating Piecewise Linear Continuous "
                   "Finite Element spatial discretizaiton.";
}
