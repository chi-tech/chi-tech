#include "pwl.h"

#include "ChiMath/UnknownManager/unknown_manager.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiTimer/chi_timer.h"


//###################################################################
/**Constructor.*/
chi_math::SpatialDiscretization_PWLD::
  SpatialDiscretization_PWLD(chi_mesh::MeshContinuumPtr& in_grid,
                             chi_math::finite_element::SetupFlags setup_flags,
                             chi_math::QuadratureOrder qorder,
                             chi_math::CoordinateSystemType in_cs_type) :
  SpatialDiscretization_PWLBase(in_grid, setup_flags, qorder,
                                SDMType::PIECEWISE_LINEAR_DISCONTINUOUS,
                                in_cs_type)
{
  chi::log.Log() << chi::program_timer.GetTimeString()
                << " Creating Piecewise Linear Discontinuous "
                   "Finite Element spatial discretizaiton.";

  chi::log.Log() << chi::program_timer.GetTimeString()
                << " Communicating partition neighbors.";

  if (setup_flags == chi_math::finite_element::COMPUTE_UNIT_INTEGRALS)
  {
    int qorder_min;
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
        throw std::invalid_argument("SpatialDiscretization_PWLD::SpatialDiscretization_PWLD : "
                                    "Unsupported coordinate system type encountered.");
    }

    if (static_cast<int>(line_quad_order_arbitrary.order) < qorder_min)
      chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLD::SpatialDiscretization_PWLD : "
        << "static_cast<int>(line_quad_order_arbitrary.order) < "
        << qorder_min << ".";

    if (static_cast<int>(tri_quad_order_arbitrary.order) < qorder_min)
      chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLD::SpatialDiscretization_PWLD : "
        << "static_cast<int>(tri_quad_order_arbitrary.order) < "
        << qorder_min << ".";

    if (static_cast<int>(quad_quad_order_arbitrary.order) < qorder_min)
      chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLD::SpatialDiscretization_PWLD : "
        << "static_cast<int>(quad_quad_order_arbitrary.order) < "
        << qorder_min << ".";

    if (static_cast<int>(tet_quad_order_arbitrary.order) < qorder_min)
      chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLD::SpatialDiscretization_PWLD : "
        << "static_cast<int>(tet_quad_order_arbitrary.order) < "
        << qorder_min << ".";
  }

  CreateCellMappings();
  if (setup_flags != chi_math::finite_element::NO_FLAGS_SET)
  {
    PreComputeCellSDValues();
    PreComputeNeighborCellSDValues();
  }
  OrderNodes();
  chi::log.Log() << chi::program_timer.GetTimeString()
                << " Done creating Piecewise Linear Discontinuous "
                   "Finite Element spatial discretizaiton.";
}


//###################################################################
/**Construct a shared object using the protected constructor.*/
std::shared_ptr<chi_math::SpatialDiscretization_PWLD>
chi_math::SpatialDiscretization_PWLD::
  New(chi_mesh::MeshContinuumPtr& in_grid,
      finite_element::SetupFlags setup_flags/*=finite_element::NO_FLAGS_SET*/,
      QuadratureOrder qorder/*=QuadratureOrder::SECOND*/,
      CoordinateSystemType in_cs_type/*=CoordinateSystemType::CARTESIAN*/)

{
  if (in_grid == nullptr)
    throw std::invalid_argument(
      "Null supplied as grid to SpatialDiscretization_PWLD.");
  return std::shared_ptr<SpatialDiscretization_PWLD>(
    new SpatialDiscretization_PWLD(in_grid, setup_flags, qorder, in_cs_type));
}
