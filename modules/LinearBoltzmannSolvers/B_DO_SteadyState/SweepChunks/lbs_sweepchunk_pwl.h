#ifndef LBS_SWEEPCHUNK_PWL_H
#define LBS_SWEEPCHUNK_PWL_H

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"

#include "B_DO_SteadyState/lbs_DO_steady_state.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

#include "ChiMesh/SweepUtilities/sweepchunk_base.h"

typedef std::map<int,std::shared_ptr<chi_physics::TransportCrossSections>> TCrossSections;

namespace lbs
{
//###################################################################
/**Sweep chunk for cartesian PWLD discretization.*/
class SweepChunkPWL : public chi_mesh::sweep_management::SweepChunk
{
protected:
  const std::shared_ptr<chi_mesh::MeshContinuum> grid_view;
  const chi_math::SpatialDiscretization& grid_fe_view;
  const std::vector<UnitCellMatrices>& unit_cell_matrices_;
  std::vector<lbs::CellLBSView>& grid_transport_view;
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
                const chi_math::SpatialDiscretization& discretization,
                const std::vector<UnitCellMatrices>& unit_cell_matrices,
                std::vector<lbs::CellLBSView>& cell_transport_views,
                std::vector<double>& destination_phi,
                std::vector<double>& destination_psi,
                const std::vector<double>& source_moments,
                LBSGroupset& in_groupset,
                const TCrossSections& in_xsections,
                int in_num_moms,
                int in_max_num_cell_dofs);

  void Sweep(chi_mesh::sweep_management::AngleSet* angle_set) override;


  struct Upwinder
  {
    chi_mesh::sweep_management::FLUDS* fluds;
    chi_mesh::sweep_management::AngleSet* angle_set;
    size_t spls_index;
    size_t angle_set_index;
    int in_face_counter;
    int preloc_face_counter;
    int out_face_counter;
    int deploc_face_counter;
    uint64_t bndry_id;
    int angle_num;
    uint64_t cell_local_id;
    int f;
    int gs_gi;
    size_t gs_ss_begin;
    bool surface_source_active;

    const double* GetUpwindPsi(int fj, bool local, bool boundary) const;
    double* GetDownwindPsi(int fi,
                           bool local,
                           bool boundary,
                           bool reflecting_bndry) const;
  };

};
}


#endif
