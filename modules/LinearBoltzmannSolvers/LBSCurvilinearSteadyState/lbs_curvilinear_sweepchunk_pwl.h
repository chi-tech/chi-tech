#ifndef LBS_CURVILINEAR_SWEEPCHUNK_PWL_H
#define LBS_CURVILINEAR_SWEEPCHUNK_PWL_H

#include "B_LBSSteadyState/SweepChunks/lbs_sweepchunk_pwl.h"
#include "B_LBSSteadyState/Groupset/lbs_groupset.h"



namespace lbs_curvilinear
{
  class SweepChunkPWL;
}


/** A sweep-chunk in point-symmetric and axial-symmetric
 *  curvilinear coordinates. */
class lbs_curvilinear::SweepChunkPWL : public lbs::SweepChunkPWL
{
//  Attributes
private:
  /** Spatial discretisation of secondary cell view (spatial discretisation
   *  of primary cell view managed by the base class). */
  chi_math::SpatialDiscretization_PWLD& grid_fe_view_secondary;
  /** Unknown manager. */
  chi_math::UnknownManager unknown_manager;
  /** Sweeping dependency angular intensity (for each polar level). */
  std::vector<double> psi_sweep;
  /** Mapping from direction linear index to direction polar level. */
  std::map<unsigned int, unsigned int> map_polar_level;
  /** Normal vector to determine symmetric boundary condition. */
  chi_mesh::Vector3 normal_vector_boundary;

//  Methods
public:
  /** Constructor. */
  SweepChunkPWL(std::shared_ptr<chi_mesh::MeshContinuum> grid_ptr,
                chi_math::SpatialDiscretization_PWLD& discretization_primary,
                chi_math::SpatialDiscretization_PWLD& discretization_secondary,
                std::vector<lbs::CellLBSView>& cell_transport_views,
                std::vector<double>& destination_phi,
                std::vector<double>& destination_psi,
                const std::vector<double>& source_moments,
                lbs::LBSGroupset& in_groupset,
                const TCrossSections& in_xsections,
                int in_num_moms,
                int in_max_num_cell_dofs);

  void Sweep(chi_mesh::sweep_management::AngleSet* angle_set) override;
};

#endif // LBS_CURVILINEAR_SWEEPCHUNK_PWL_H
