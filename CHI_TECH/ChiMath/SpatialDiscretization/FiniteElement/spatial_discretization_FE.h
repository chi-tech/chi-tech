#ifndef SPATIAL_DISCRETIZATION_FE_H
#define SPATIAL_DISCRETIZATION_FE_H

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/UnknownManager/unknown_manager.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

//###################################################################
/**Base Finite Element spatial discretization class.
 * */
class SpatialDiscretization_FE : public SpatialDiscretization
{
protected:
  std::vector<chi_math::finite_element::UnitIntegralData> fe_unit_integrals;
  std::vector<chi_math::finite_element::InternalQuadraturePointData> fe_vol_qp_data;
  std::vector<std::vector<chi_math::finite_element::FaceQuadraturePointData>> fe_srf_qp_data;

  bool integral_data_initialized=false;
  bool qp_data_initialized=false;

protected:
  SpatialDiscretization_FE(int dim,
                           chi_mesh::MeshContinuumPtr in_grid,
                           SDMType in_type =
                           SDMType::UNDEFINED) :
    SpatialDiscretization(dim,in_grid,in_type)
  {}

public:
  virtual
  const chi_math::finite_element::UnitIntegralData&
  GetUnitIntegrals(chi_mesh::Cell& cell) const
  {
    if (not integral_data_initialized)
      throw std::invalid_argument("SpatialDiscretization_FE::GetUnitIntegrals "
                                  "called without integrals being initialized.");
    return fe_unit_integrals[cell.local_id];
  }

  virtual
  const chi_math::finite_element::InternalQuadraturePointData&
  GetQPData_Volumetric(chi_mesh::Cell& cell) const
  {
    if (not qp_data_initialized)
      throw std::invalid_argument("SpatialDiscretization_FE::GetUnitIntegrals "
                                  "called without integrals being initialized.");
    return fe_vol_qp_data[cell.local_id];
  }

  virtual
  const chi_math::finite_element::FaceQuadraturePointData&
  GetQPData_Surface(const chi_mesh::Cell& cell,
                       const unsigned int face) const
  {
    if (not qp_data_initialized)
      throw std::invalid_argument("SpatialDiscretization_FE::GetUnitIntegrals "
                                  "called without integrals being initialized.");
    return fe_srf_qp_data[cell.local_id][face];
  }

  virtual ~SpatialDiscretization_FE() = default;
};


#endif