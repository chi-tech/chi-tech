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

    std::vector<UIData>                  fe_unit_integrals_;
    std::vector<QPDataVol>               fe_vol_qp_data_;
    std::vector<std::vector<QPDataFace>> fe_srf_qp_data_;

    std::map<uint64_t, UIData>                  nb_fe_unit_integrals_;
    std::map<uint64_t, QPDataVol>               nb_fe_vol_qp_data_;
    std::map<uint64_t, std::vector<QPDataFace>> nb_fe_srf_qp_data_;

    bool integral_data_initialized_ = false;
    bool qp_data_initialized_ = false;

    bool nb_integral_data_initialized_ = false;
    bool nb_qp_data_initialized_ = false;

    UIData     scratch_intgl_data_;
    QPDataVol  scratch_vol_qp_data_;
    QPDataFace scratch_face_qp_data_;

    const finite_element::SetupFlags setup_flags_;


  protected:
    explicit
    SpatialDiscretization_FE(const chi_mesh::MeshContinuum& in_grid,
                             CoordinateSystemType in_cs_type =
                             CoordinateSystemType::CARTESIAN,
                             SDMType in_type =
                             SDMType::UNDEFINED,
                             finite_element::SetupFlags in_setup_flags=
                             finite_element::NO_FLAGS_SET) :
        SpatialDiscretization(in_grid, in_cs_type, in_type),
        setup_flags_(in_setup_flags)
    {}

  public:
    virtual
    const finite_element::UnitIntegralData&
      GetUnitIntegrals(const chi_mesh::Cell& cell)
    {
      if (not integral_data_initialized_)
        throw std::invalid_argument("SpatialDiscretization_FE::GetUnitIntegrals "
                                    "called without integrals being initialized."
                                    " Set flag COMPUTE_UNIT_INTEGRALS.");
      return fe_unit_integrals_[cell.local_id_];
    }

    virtual
    const finite_element::InternalQuadraturePointData&
      GetQPData_Volumetric(const chi_mesh::Cell& cell)
    {
      if (not qp_data_initialized_)
        throw std::invalid_argument("SpatialDiscretization_FE::GetQPData_Volumetric "
                                    "called without integrals being initialized."
                                    " Set flag INIT_QP_DATA.");
      return fe_vol_qp_data_[cell.local_id_];
    }

    virtual
    const finite_element::FaceQuadraturePointData&
      GetQPData_Surface(const chi_mesh::Cell& cell,
                        const unsigned int face)
    {
      if (not qp_data_initialized_)
        throw std::invalid_argument("SpatialDiscretization_FE::GetQPData_Surface "
                                    "called without integrals being initialized."
                                    " Set flag INIT_QP_DATA.");
      return fe_srf_qp_data_[cell.local_id_][face];
    }

    virtual ~SpatialDiscretization_FE() = default;
  };
}


#endif