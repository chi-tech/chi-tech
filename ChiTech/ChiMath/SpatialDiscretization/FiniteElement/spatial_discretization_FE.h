#ifndef SPATIAL_DISCRETIZATION_FE_H
#define SPATIAL_DISCRETIZATION_FE_H

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/UnknownManager/unknown_manager.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

//###################################################################
namespace chi_math
{
  /**Base Finite Element spatial discretization class.
     * */
  class SpatialDiscretization_FE : public chi_math::SpatialDiscretization
  {
  protected:
    typedef finite_element::UnitIntegralData UIData;
    typedef finite_element::InternalQuadraturePointData QPDataVol;
    typedef finite_element::FaceQuadraturePointData QPDataFace;

    std::vector<UIData>                  fe_unit_integrals;
    std::vector<QPDataVol>               fe_vol_qp_data;
    std::vector<std::vector<QPDataFace>> fe_srf_qp_data;

    bool integral_data_initialized=false;
    bool qp_data_initialized=false;

    const finite_element::SetupFlags setup_flags;

  protected:
    SpatialDiscretization_FE(int dim,
                             chi_mesh::MeshContinuumPtr& in_grid,
                             CoordinateSystemType in_cs_type =
                             CoordinateSystemType::CARTESIAN,
                             SDMType in_type =
                             SDMType::UNDEFINED,
                             finite_element::SetupFlags in_setup_flags=
                             finite_element::NO_FLAGS_SET) :
      SpatialDiscretization(dim, in_grid, in_cs_type, in_type),
      setup_flags(in_setup_flags)
    {}

  public:
    virtual
    const finite_element::UnitIntegralData&
      GetUnitIntegrals(const chi_mesh::Cell& cell)
    {
      if (not integral_data_initialized)
        throw std::invalid_argument("SpatialDiscretization_FE::GetUnitIntegrals "
                                    "called without integrals being initialized."
                                    " Set flag COMPUTE_UNIT_INTEGRALS.");
      return fe_unit_integrals[cell.local_id];
    }

    virtual
    const finite_element::InternalQuadraturePointData&
      GetQPData_Volumetric(const chi_mesh::Cell& cell)
    {
      if (not qp_data_initialized)
        throw std::invalid_argument("SpatialDiscretization_FE::GetQPData_Volumetric "
                                    "called without integrals being initialized."
                                    " Set flag INIT_QP_DATA.");
      return fe_vol_qp_data[cell.local_id];
    }

    virtual
    const finite_element::FaceQuadraturePointData&
      GetQPData_Surface(const chi_mesh::Cell& cell,
                        const unsigned int face)
    {
      if (not qp_data_initialized)
        throw std::invalid_argument("SpatialDiscretization_FE::GetQPData_Surface "
                                    "called without integrals being initialized."
                                    " Set flag INIT_QP_DATA.");
      return fe_srf_qp_data[cell.local_id][face];
    }

    virtual ~SpatialDiscretization_FE() = default;
  };
}


#endif