#ifndef CHITECH_LBS_SWEEPCHUNK_H
#define CHITECH_LBS_SWEEPCHUNK_H

#include "ChiMesh/SweepUtilities/sweepchunk_base.h"

#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

namespace lbs
{

// ##################################################################
/**Simple utility structure for controll counters and calls
 * to upstream data.*/
struct SweepSurfaceStatusInfo
{
  chi_mesh::sweep_management::AngleSet* angle_set;
  chi_mesh::sweep_management::FLUDS* fluds;
  bool surface_source_active = false;
  size_t gs_ss_size_ = 0;
  size_t gs_ss_begin_ = 0;
  int gs_gi_ = 0;

  size_t spls_index = 0;
  uint64_t cell_local_id = 0;

  size_t angle_set_index = 0;
  size_t angle_num = 0;

  int in_face_counter = 0;
  int preloc_face_counter = 0;
  int out_face_counter = 0;
  int deploc_face_counter = 0;

  uint64_t bndry_id = 0;
  int f = 0;

  bool on_local_face = false;
  bool on_boundary = false;
  bool is_reflecting_bndry_ = false;

  const double* GetUpwindPsi(int fj) const;
  double* GetDownwindPsi(int fi) const;
};

// ##################################################################
/**The new sweep chunk class.*/
class LBSSweepChunk : public chi_mesh::sweep_management::SweepChunk
{
protected:
  typedef std::function<void()> CallbackFunction;

protected:
  const chi_mesh::MeshContinuum& grid_;
  const chi_math::SpatialDiscretization& grid_fe_view_;
  const std::vector<UnitCellMatrices>& unit_cell_matrices_;
  std::vector<lbs::CellLBSView>& grid_transport_view_;
  const std::vector<double>& q_moments_;
  LBSGroupset& groupset_;
  const std::map<int, XSPtr>& xs_;
  const int num_moments_;
  const size_t num_groups_;
  const bool save_angular_flux_;

  size_t gs_ss_size_ = 0;
  size_t gs_ss_begin_ = 0;
  int gs_gi_ = 0;

  // Runtime params
  std::vector<std::vector<double>> Amat_;
  std::vector<std::vector<double>> Atemp_;
  std::vector<double> source_;
  std::vector<std::vector<double>> b_;

  SweepSurfaceStatusInfo sweep_surface_status_info_;

  uint64_t cell_local_id_ = 0;
  const chi_mesh::Cell* cell_ = nullptr;
  const chi_math::CellMapping* cell_mapping_ = nullptr;
  CellLBSView* cell_transport_view_ = nullptr;
  size_t cell_num_faces_ = 0;
  size_t cell_num_nodes_ = 0;
  const MatVec3* G_ = nullptr;
  const MatDbl* M_ = nullptr;
  const std::vector<MatDbl>* M_surf_ = nullptr;
  const std::vector<VecDbl>* IntS_shapeI_ = nullptr;

  /**Callbacks at phase 1 : cell data established*/
  std::vector<CallbackFunction> cell_data_callbacks_;

  std::vector<bool> face_incident_flags_;
  std::vector<double> face_mu_values_;
  size_t direction_num_ = 0;
  chi_mesh::Vector3 omega_;
  double direction_qweight_ = 0.0;

  /**Callbacks at phase 2 : direction data established*/
  std::vector<CallbackFunction> direction_data_callbacks_and_kernels_;

  /**Callbacks at phase 3 : Surface integrals*/
  std::vector<CallbackFunction> surface_integral_kernels_;

  size_t g_ = 0;
  size_t gsg_ = 0;
  double sigma_tg_ = 0.0;

  /**Callbacks at phase 4 : group by group mass terms*/
  std::vector<CallbackFunction> mass_term_kernels_;

  /**Callbacks at phase 5 : flux updates*/
  std::vector<CallbackFunction> flux_update_kernels_;

  /**Callbacks at phase 6 : Post cell-dir sweep*/
  std::vector<CallbackFunction> post_cell_dir_sweep_callbacks_;

private:
  std::map<std::string, CallbackFunction> kernels_;

public:
  LBSSweepChunk(const chi_mesh::MeshContinuum& grid,
                const chi_math::SpatialDiscretization& discretization,
                const std::vector<UnitCellMatrices>& unit_cell_matrices,
                std::vector<lbs::CellLBSView>& cell_transport_views,
                std::vector<double>& destination_phi,
                std::vector<double>& destination_psi,
                const std::vector<double>& source_moments,
                LBSGroupset& groupset,
                const std::map<int, XSPtr>& xs,
                int num_moments,
                int max_num_cell_dofs);

  // 01
  void Sweep(chi_mesh::sweep_management::AngleSet* angle_set) override;

protected:
  // 02 operations
  void RegisterKernel(const std::string& name, CallbackFunction function);
  CallbackFunction Kernel(const std::string& name) const;
  static void ExecuteKernels(const std::vector<CallbackFunction>& kernels);
  virtual void OutgoingSurfaceOperations();

  // 03 kernels
  void KernelFEMVolumetricGradientTerm();
  void KernelFEMUpwindSurfaceIntegrals();
  void KernelFEMSTDMassTerms();
  void KernelPhiUpdate();
  void KernelPsiUpdate();
};

} // namespace lbs

#endif // CHITECH_LBS_SWEEPCHUNK_H
