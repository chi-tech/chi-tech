#include <utility>

#include "lbs_sweepchunk.h"

#include "chi_log_exceptions.h"

namespace lbs
{

// ##################################################################
void LBSSweepChunk::RegisterKernel(const std::string& name,
                                   CallbackFunction function)
{
  ChiInvalidArgumentIf(kernels_.count(name) > 0,
                       "Attempting to register kernel with name \"" + name +
                         "\" but the kernel already exists.");

  kernels_[name] = std::move(function);
}

// ##################################################################
/**Returns a kernel if the given name exists.*/
LBSSweepChunk::CallbackFunction
LBSSweepChunk::Kernel(const std::string& name) const
{
  ChiInvalidArgumentIf(kernels_.count(name) == 0,
                       "No register kernel with name \"" + name + "\" found");
  return kernels_.at(name);
}

// ##################################################################
/**Executes the supplied kernels list.*/
void LBSSweepChunk::ExecuteKernels(const std::vector<CallbackFunction>& kernels)
{
  for (auto& kernel : kernels)
    kernel();
}

// ##################################################################
/**Operations when outgoing fluxes are handled including passing
 * face angular fluxes downstream and computing
 * balance parameters (i.e. outflow)
 * */
void LBSSweepChunk::OutgoingSurfaceOperations()
{
  const size_t f = sweep_surface_status_info_.f;
  const auto& IntF_shapeI = (*IntS_shapeI_)[f];
  const double mu = face_mu_values_[f];
  const double wt = direction_qweight_;

  const auto& sss_info = sweep_surface_status_info_;

  const size_t num_face_nodes = cell_mapping_->NumFaceNodes(f);
  for (int fi = 0; fi < num_face_nodes; ++fi)
  {
    const int i = cell_mapping_->MapFaceNode(f, fi);

    double* psi = sweep_surface_status_info_.GetDownwindPsi(fi);

    if (not sss_info.on_boundary or sss_info.is_reflecting_bndry_)
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        psi[gsg] = b_[gsg][i];
    if (sss_info.on_boundary and not sss_info.is_reflecting_bndry_)
      for (int gsg = 0; gsg < gs_ss_size_; ++gsg)
        cell_transport_view_->AddOutflow(gs_gi_ + gsg,
                                         wt * mu * b_[gsg][i] * IntF_shapeI[i]);
  } // for fi
}

} // namespace lbs