#ifndef SPATIAL_DISCRETIZATION_FE_H
#define SPATIAL_DISCRETIZATION_FE_H

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

//###################################################################
/**Base Finite Element spatial discretization class.
 * */
class SpatialDiscretization_FE : public SpatialDiscretization
{
protected:
  SpatialDiscretization_FE(int dim,
                           chi_mesh::MeshContinuumPtr in_grid,
                           SDMType in_type =
                           SDMType::UNDEFINED) :
    SpatialDiscretization(dim,in_grid,in_type)
  {}

public:
  virtual ~SpatialDiscretization_FE() = default;
};


#endif