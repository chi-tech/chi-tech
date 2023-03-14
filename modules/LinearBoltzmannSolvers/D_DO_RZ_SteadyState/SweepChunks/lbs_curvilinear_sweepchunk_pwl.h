#ifndef LBS_CURVILINEAR_SWEEPCHUNK_PWL_H
#define LBS_CURVILINEAR_SWEEPCHUNK_PWL_H

#include "B_DO_Solver/SweepChunks/lbs_sweepchunk_pwl.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"



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
  const std::vector<lbs::UnitCellMatrices>& secondary_unit_cell_matrices_;
  /** Unknown manager. */
  chi_math::UnknownManager unknown_manager_;
  /** Sweeping dependency angular intensity (for each polar level). */
  std::vector<double> psi_sweep_;
  /** Mapping from direction linear index to direction polar level. */
  std::map<unsigned int, unsigned int> map_polar_level_;
  /** Normal vector to determine symmetric boundary condition. */
  chi_mesh::Vector3 normal_vector_boundary_;

//  Methods
public:
  /** Constructor. */
  SweepChunkPWL(std::shared_ptr<chi_mesh::MeshContinuum> grid_ptr,
                chi_math::SpatialDiscretization& discretization_primary,
                const std::vector<lbs::UnitCellMatrices>& unit_cell_matrices,
                const std::vector<lbs::UnitCellMatrices>& secondary_unit_cell_matrices,
                std::vector<lbs::CellLBSView>& cell_transport_views,
                std::vector<double>& destination_phi,
                std::vector<double>& destination_psi,
                const std::vector<double>& source_moments,
                lbs::LBSGroupset& groupset,
                const std::map<int, lbs::XSPtr>& xs,
                int num_moments,
                int max_num_cell_dofs);

  void Sweep(chi_mesh::sweep_management::AngleSet* angle_set) override;
};

#endif // LBS_CURVILINEAR_SWEEPCHUNK_PWL_H
