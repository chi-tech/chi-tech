#ifndef CHITECH_LAGRANGEBASE_H
#define CHITECH_LAGRANGEBASE_H

#include "math/SpatialDiscretization/FiniteElement/FiniteElementBase.h"

namespace chi_math::spatial_discretization
{

/**Base class for Lagrange spatial discretizations.
* \ingroup doc_SpatialDiscretization*/
class LagrangeBase : public FiniteElementBase
{
protected:
  LagrangeBase(const chi_mesh::MeshContinuum& grid,
               QuadratureOrder q_order,
               SDMType sdm_type,
               CoordinateSystemType cs_type);

  void CreateCellMappings();
};

}

#endif // CHITECH_LAGRANGEBASE_H
