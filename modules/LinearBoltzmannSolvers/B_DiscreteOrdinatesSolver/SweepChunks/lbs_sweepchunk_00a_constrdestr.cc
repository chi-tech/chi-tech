#include "lbs_sweepchunk.h"

namespace lbs
{

LBSSweepChunk::LBSSweepChunk(
  const chi_mesh::MeshContinuum& grid,
  const chi_math::SpatialDiscretization& discretization,
  const std::vector<UnitCellMatrices>& unit_cell_matrices,
  std::vector<lbs::CellLBSView>& cell_transport_views,
  std::vector<double>& destination_phi,
  std::vector<double>& destination_psi,
  const std::vector<double>& source_moments,
  const LBSGroupset& groupset,
  const std::map<int, XSPtr>& xs,
  int num_moments,
  int max_num_cell_dofs)
  : SweepChunk(destination_phi, destination_psi),
    grid_(grid),
    grid_fe_view_(discretization),
    unit_cell_matrices_(unit_cell_matrices),
    grid_transport_view_(cell_transport_views),
    q_moments_(source_moments),
    groupset_(groupset),
    xs_(xs),
    num_moments_(num_moments),
    num_groups_(groupset.groups_.size()),
    save_angular_flux_(!destination_psi.empty())
{
  Amat_.resize(max_num_cell_dofs, std::vector<double>(max_num_cell_dofs));
  Atemp_.resize(max_num_cell_dofs, std::vector<double>(max_num_cell_dofs));
  b_.resize(num_groups_, std::vector<double>(max_num_cell_dofs, 0.0));
  source_.resize(max_num_cell_dofs, 0.0);

  // ================================== Register kernels
  RegisterKernel("FEMVolumetricGradTerm",
    std::bind(&LBSSweepChunk::KernelFEMVolumetricGradientTerm, this));
  RegisterKernel("FEMUpwindSurfaceIntegrals",
    std::bind(&LBSSweepChunk::KernelFEMUpwindSurfaceIntegrals, this));
  RegisterKernel("FEMSSTDMassTerms",
    std::bind(&LBSSweepChunk::KernelFEMSTDMassTerms, this));
  RegisterKernel("KernelPhiUpdate",
    std::bind(&LBSSweepChunk::KernelPhiUpdate, this));
  RegisterKernel("KernelPsiUpdate",
    std::bind(&LBSSweepChunk::KernelPsiUpdate, this));

  // ================================== Setup callbacks
  cell_data_callbacks_ = {};

  direction_data_callbacks_and_kernels_ = {
    Kernel("FEMVolumetricGradTerm")};

  surface_integral_kernels_ = {Kernel("FEMUpwindSurfaceIntegrals")};

  mass_term_kernels_ = {Kernel("FEMSSTDMassTerms")};

  flux_update_kernels_ = {Kernel("KernelPhiUpdate"),
                          Kernel("KernelPsiUpdate")};

  post_cell_dir_sweep_callbacks_ = {};
}

} // namespace lbs