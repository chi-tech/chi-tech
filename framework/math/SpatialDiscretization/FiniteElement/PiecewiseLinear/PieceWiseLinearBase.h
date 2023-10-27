#ifndef SPATIAL_DISCRETIZATION_PWL_BASE_H
#define SPATIAL_DISCRETIZATION_PWL_BASE_H

#include "math/SpatialDiscretization/FiniteElement/FiniteElementBase.h"

namespace chi_math::spatial_discretization
{

/**Base class for PieceWiseLinear based discretization.
* \ingroup doc_SpatialDiscretization*/
class PieceWiseLinearBase : public FiniteElementBase
{
protected:
  explicit PieceWiseLinearBase(const chi_mesh::MeshContinuum& grid,
                               QuadratureOrder q_order,
                               SDMType sdm_type,
                               CoordinateSystemType cs_type);

  void CreateCellMappings();
};

} // namespace chi_math::spatial_discretization

#endif // SPATIAL_DISCRETIZATION_PWL_BASE_H
