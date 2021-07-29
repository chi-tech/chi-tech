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

  int qorder_min = static_cast<int>(chi_math::QuadratureOrder::INVALID_ORDER);
  switch (cs_type)
  {
    case chi_math::CoordinateSystemType::CARTESIAN:
    {
      qorder_min = static_cast<int>(chi_math::QuadratureOrder::SECOND);
      break;
    }
    case chi_math::CoordinateSystemType::CYLINDRICAL:
    {
      qorder_min = static_cast<int>(chi_math::QuadratureOrder::THIRD);
      break;
    }
    case chi_math::CoordinateSystemType::SPHERICAL:
    {
      qorder_min = static_cast<int>(chi_math::QuadratureOrder::FOURTH);
      break;
    }
    default:
      throw std::invalid_argument("SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
                                  "Unsupported coordinate system type encountered.");
  }

  if (static_cast<int>(line_quad_order_arbitrary.order) < qorder_min)
    chi_log.Log(LOG_ALLWARNING)
      << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
      << "static_cast<int>(line_quad_order_arbitrary.order) < "
      << qorder_min << ".";

  if (static_cast<int>(tri_quad_order_arbitrary.order) < qorder_min)
    chi_log.Log(LOG_ALLWARNING)
      << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
      << "static_cast<int>(tri_quad_order_arbitrary.order) < "
      << qorder_min << ".";

  if (static_cast<int>(quad_quad_order_arbitrary.order) < qorder_min)
    chi_log.Log(LOG_ALLWARNING)
      << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
      << "static_cast<int>(quad_quad_order_arbitrary.order) < "
      << qorder_min << ".";

  if (static_cast<int>(tet_quad_order_arbitrary.order) < qorder_min)
    chi_log.Log(LOG_ALLWARNING)
      << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
      << "static_cast<int>(tet_quad_order_arbitrary.order) < "
      << qorder_min << ".";

  if (static_cast<int>(hex_quad_order_arbitrary.order) < qorder_min)
    chi_log.Log(LOG_ALLWARNING)
      << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
      << "static_cast<int>(hex_quad_order_arbitrary.order) < "
      << qorder_min << ".";

  PreComputeCellSDValues();

  OrderNodes();
  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Done creating Piecewise Linear Continuous "
                   "Finite Element spatial discretizaiton.";
}
