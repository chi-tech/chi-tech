#ifndef LBS_CURVILINEAR_SWEEPCHUNK_PWL_H
#define LBS_CURVILINEAR_SWEEPCHUNK_PWL_H

#include "LinearBoltzmannSolver/SweepChunks/lbs_sweepchunk_pwl.h"


namespace LBSCurvilinear
{
  class SweepChunkPWL;
}


/** A sweep-chunk in point-symmetric and axial-symmetric
 *  curvilinear coordinates. */
class LBSCurvilinear::SweepChunkPWL : public LinearBoltzmann::LBSSweepChunkPWL
{
//  Attributes
private:
  /** Spatial discretisation of secondary cell view (spatial discretisation
   *  of primary cell view managed by the base class). */
  SpatialDiscretization_PWLD& grid_fe_view_secondary;
  /** Unknown manager. */
  chi_math::UnknownManager unknown_manager;
  /** Starting direction angular intensity (for each polar level). */
  std::vector<double> psi_start;
  /** Sweeping dependency angular intensity (for each polar level). */
  std::vector<double> psi_sweep;
  /** Mapping from direction linear index to direction polar level. */
  std::map<unsigned int, unsigned int> map_polar_level;
  /** Mapping from direction linear index to information on start/final. */
  std::map<unsigned int, std::pair<bool, bool>> map_start_final;
  /** Normal vector to determine symmetric boundary condition. */
  chi_mesh::Vector3 normal_vector_boundary;

//  Methods
public:
  /** Constructor. */
  SweepChunkPWL(std::shared_ptr<chi_mesh::MeshContinuum> grid_ptr,
                SpatialDiscretization_PWLD& discretization_primary,
                SpatialDiscretization_PWLD& discretization_secondary,
                std::vector<LinearBoltzmann::CellLBSView>& cell_transport_views,
                std::vector<double>& destination_phi,
                const std::vector<double>& source_moments,
                LBSGroupset& in_groupset,
                const TCrossSections& in_xsections,
                const int in_num_moms,
                const int in_max_num_cell_dofs);

  void Sweep(chi_mesh::sweep_management::AngleSet* angle_set) override;
};

#endif // LBS_CURVILINEAR_SWEEPCHUNK_PWL_H
