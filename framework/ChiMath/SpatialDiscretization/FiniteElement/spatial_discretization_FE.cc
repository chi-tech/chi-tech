#include "spatial_discretization_FE.h"

chi_math::finite_element::SetupFlags
chi_math::SpatialDiscretization_FE::GetSetupFlags() const
{
  return setup_flags_;
}

chi_math::QuadratureOrder
  chi_math::SpatialDiscretization_FE::GetQuadratureOrder() const
{
  return q_order_;
}