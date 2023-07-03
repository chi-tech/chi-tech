#include "fv.h"

#include "chi_runtime.h"
#include "chi_log_exceptions.h"

// ##################################################################
/**Only constructor for this method.*/
chi_math::SpatialDiscretization_FV::SpatialDiscretization_FV(
  const chi_mesh::MeshContinuum& in_grid,
  chi_math::CoordinateSystemType in_cs_type)
  : chi_math::SpatialDiscretization(in_grid, in_cs_type, SDMType::FINITE_VOLUME)
{
  CreateCellMappings();

  OrderNodes();
}

// ##################################################################
/**Publicly accessible construction handler.*/
std::shared_ptr<chi_math::SpatialDiscretization_FV>
chi_math::SpatialDiscretization_FV::New(const chi_mesh::MeshContinuum& in_grid,
    chi_math::CoordinateSystemType in_cs_type/* =
      chi_math::CoordinateSystemType::CARTESIAN*/)
{
  // First try to find an existing spatial discretization that matches the
  // one requested.
  for (auto& sdm : Chi::sdm_stack)
    if (sdm->Type() == SpatialDiscretizationType::FINITE_VOLUME and
        std::addressof(sdm->Grid()) == std::addressof(in_grid) and
        sdm->GetCoordinateSystemType() == in_cs_type)
    {
      auto sdm_ptr = std::dynamic_pointer_cast<SpatialDiscretization_FV>(sdm);

      ChiLogicalErrorIf(not sdm_ptr, "Casting failure");

      return sdm_ptr;
    }

  // If no existing discretization was found then go ahead and make a
  // new one
  auto new_sdm = std::shared_ptr<chi_math::SpatialDiscretization_FV>(
    new SpatialDiscretization_FV(in_grid, in_cs_type));

  Chi::sdm_stack.push_back(new_sdm);

  return new_sdm;
}