#include "pwlc.h"

#include "chi_log.h"

#include "ChiTimer/chi_timer.h"




//###################################################################
/**Constructor.*/
chi_math::SpatialDiscretization_PWLC::
  SpatialDiscretization_PWLC(chi_mesh::MeshContinuumPtr& in_grid,
                             chi_math::finite_element::SetupFlags setup_flags,
                             chi_math::QuadratureOrder qorder,
                             chi_math::CoordinateSystemType in_cs_type) :
  SpatialDiscretization_FE(in_grid, in_cs_type,
                           SDMType::PIECEWISE_LINEAR_CONTINUOUS,
                           setup_flags),
  line_quad_order_arbitrary(qorder),
  tri_quad_order_arbitrary(qorder),
  quad_quad_order_arbitrary(qorder),
  tet_quad_order_arbitrary(qorder),
  hex_quad_order_arbitrary(qorder)
{
  chi::log.Log() << chi::program_timer.GetTimeString()
                << " Creating Piecewise Linear Continuous "
                   "Finite Element spatial discretizaiton.";

  if (setup_flags == chi_math::finite_element::COMPUTE_UNIT_INTEGRALS)
  {
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
        throw std::invalid_argument(
          "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
          "Unsupported coordinate system type encountered.");
    }

    if (static_cast<int>(line_quad_order_arbitrary.order) < qorder_min)
      chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
        << "static_cast<int>(line_quad_order_arbitrary.order) < "
        << qorder_min << ".";

    if (static_cast<int>(tri_quad_order_arbitrary.order) < qorder_min)
      chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
        << "static_cast<int>(tri_quad_order_arbitrary.order) < "
        << qorder_min << ".";

    if (static_cast<int>(quad_quad_order_arbitrary.order) < qorder_min)
      chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
        << "static_cast<int>(quad_quad_order_arbitrary.order) < "
        << qorder_min << ".";

    if (static_cast<int>(tet_quad_order_arbitrary.order) < qorder_min)
      chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
        << "static_cast<int>(tet_quad_order_arbitrary.order) < "
        << qorder_min << ".";

    if (static_cast<int>(hex_quad_order_arbitrary.order) < qorder_min)
      chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLC::SpatialDiscretization_PWLC : "
        << "static_cast<int>(hex_quad_order_arbitrary.order) < "
        << qorder_min << ".";
  }

  if (setup_flags != chi_math::finite_element::NO_FLAGS_SET)
    PreComputeCellSDValues();

  CreateCellMappings();
  OrderNodes();
  chi::log.Log() << chi::program_timer.GetTimeString()
                << " Done creating Piecewise Linear Continuous "
                   "Finite Element spatial discretizaiton.";
}
