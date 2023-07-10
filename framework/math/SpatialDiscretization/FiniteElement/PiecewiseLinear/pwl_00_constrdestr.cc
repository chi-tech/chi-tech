#include "pwl.h"

#include "math/UnknownManager/unknown_manager.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "utils/chi_timer.h"

// ###################################################################
/**Constructor.*/
chi_math::SpatialDiscretization_PWLD::SpatialDiscretization_PWLD(
  const chi_mesh::MeshContinuum& in_grid,
  chi_math::finite_element::SetupFlags setup_flags,
  chi_math::QuadratureOrder qorder,
  chi_math::CoordinateSystemType in_cs_type)
  : SpatialDiscretization_PWLBase(in_grid,
                                  setup_flags,
                                  qorder,
                                  SDMType::PIECEWISE_LINEAR_DISCONTINUOUS,
                                  in_cs_type)
{
  if (setup_flags == chi_math::finite_element::COMPUTE_UNIT_INTEGRALS)
  {
    int qorder_min;
    switch (coord_sys_type_)
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
          "SpatialDiscretization_PWLD::SpatialDiscretization_PWLD : "
          "Unsupported coordinate system type encountered.");
    }

    if (static_cast<int>(line_quad_order_arbitrary_.order_) < qorder_min)
      Chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLD::SpatialDiscretization_PWLD : "
        << "static_cast<int>(line_quad_order_arbitrary.order) < " << qorder_min
        << ".";

    if (static_cast<int>(tri_quad_order_arbitrary_.order_) < qorder_min)
      Chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLD::SpatialDiscretization_PWLD : "
        << "static_cast<int>(tri_quad_order_arbitrary.order) < " << qorder_min
        << ".";

    if (static_cast<int>(quad_quad_order_arbitrary_.order_) < qorder_min)
      Chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLD::SpatialDiscretization_PWLD : "
        << "static_cast<int>(quad_quad_order_arbitrary.order) < " << qorder_min
        << ".";

    if (static_cast<int>(tet_quad_order_arbitrary_.order_) < qorder_min)
      Chi::log.LogAllWarning()
        << "SpatialDiscretization_PWLD::SpatialDiscretization_PWLD : "
        << "static_cast<int>(tet_quad_order_arbitrary.order) < " << qorder_min
        << ".";
  }

  CreateCellMappings();
  if (setup_flags != chi_math::finite_element::NO_FLAGS_SET)
  {
    PreComputeCellSDValues();
    PreComputeNeighborCellSDValues();
  }
  OrderNodes();
}

// ###################################################################
/**Construct a shared object using the protected constructor.*/
std::shared_ptr<chi_math::SpatialDiscretization_PWLD>
chi_math::SpatialDiscretization_PWLD::New(
  const chi_mesh::MeshContinuum& in_grid,
  finite_element::SetupFlags setup_flags /*=finite_element::NO_FLAGS_SET*/,
  QuadratureOrder qorder /*=QuadratureOrder::SECOND*/,
  CoordinateSystemType in_cs_type /*=CoordinateSystemType::CARTESIAN*/)

{
  const auto PWLD = SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS;
  // First try to find an existing spatial discretization that matches the
  // one requested.
  for (auto& sdm : Chi::sdm_stack)
    if (sdm->Type() == PWLD and
        std::addressof(sdm->Grid()) == std::addressof(in_grid) and
        sdm->GetCoordinateSystemType() == in_cs_type)
    {
      auto fe_ptr = std::dynamic_pointer_cast<SpatialDiscretization_FE>(sdm);

      ChiLogicalErrorIf(not fe_ptr, "Casting failure to FE");

      if (fe_ptr->GetSetupFlags() != setup_flags) break;
      if (fe_ptr->GetQuadratureOrder() != qorder) break;

      auto sdm_ptr =
        std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(fe_ptr);

      ChiLogicalErrorIf(not sdm_ptr, "Casting failure");

      return sdm_ptr;
    }

  auto new_sdm = std::shared_ptr<SpatialDiscretization_PWLD>(
    new SpatialDiscretization_PWLD(in_grid, setup_flags, qorder, in_cs_type));

  Chi::sdm_stack.push_back(new_sdm);

  return new_sdm;
}
