#ifndef LBS_SWEEPCHUNK_PWL_H
#define LBS_SWEEPCHUNK_PWL_H

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiPhysics/chi_physics.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include "LinearBoltzmannSolver/lbs_linear_boltzmann_solver.h"

typedef std::vector<std::shared_ptr<chi_physics::TransportCrossSections>> TCrossSections;

namespace LinearBoltzmann
{
//###################################################################
/**Sweep chunk for cartesian PWLD discretization.*/
class SweepChunkPWL : public chi_mesh::sweep_management::SweepChunk
{
protected:
  const std::shared_ptr<chi_mesh::MeshContinuum> grid_view;
  SpatialDiscretization_PWLD& grid_fe_view;
  std::vector<LinearBoltzmann::CellLBSView>& grid_transport_view;
  std::vector<double>& psi_new_local;
  const std::vector<double>& q_moments;
  LBSGroupset& groupset;
  const TCrossSections& xsections;
  const int num_moms;
  const size_t num_grps;
  const int max_num_cell_dofs;
  const bool save_angular_flux;

  //Runtime params
  bool a_and_b_initialized;
  std::vector<std::vector<double>> Amat;
  std::vector<std::vector<double>> Atemp;
  std::vector<double> source;

public:
  std::vector<std::vector<double>> b;

  SweepChunkPWL(std::shared_ptr<chi_mesh::MeshContinuum> grid_ptr,
                SpatialDiscretization_PWLD& discretization,
                std::vector<LinearBoltzmann::CellLBSView>& cell_transport_views,
                std::vector<double>& destination_phi,
                std::vector<double>& destination_psi,
                const std::vector<double>& source_moments,
                LBSGroupset& in_groupset,
                const TCrossSections& in_xsections,
                int in_num_moms,
                int in_max_num_cell_dofs);

  void Sweep(chi_mesh::sweep_management::AngleSet* angle_set) override;

  void ZeroDestinationPsi() override
    {psi_new_local.assign(psi_new_local.size(), 0.0);}

  void ZeroIncomingDelayedPsi() override
    {groupset.angle_agg.ZeroIncomingDelayedPsi();}

  void ZeroOutgoingDelayedPsi() override
    {groupset.angle_agg.ZeroOutgoingDelayedPsi();}
};
}


#endif
