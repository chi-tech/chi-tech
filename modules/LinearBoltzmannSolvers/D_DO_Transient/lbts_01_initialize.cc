#include "lbts_transient_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Transient solver initialize routine.*/
void lbs::DiscOrdTransientSolver::Initialize()
{
  chi::log.Log() << "Initializing " << TextName() << ".";
  options_.save_angular_flux = true;
  DiscOrdKEigenvalueSolver::Initialize();
  DiscOrdKEigenvalueSolver::Execute();

  //======================================== Scale fission data
  //TODO: Determine a better methodology to handle fission scaling.

  // NOTE: This is done to ensure consistency between cross sections
  //       that may be swapped mid-simulation. For example, if one
  //       seeks to swap to cross sections with more or less absorption,
  //       then if this loop is over material_xs instead of the global
  //       stack, the two cross section sets will have different fission
  //       cross sections despite that not being intended.
  // NOTE: A potentially better way to handle this is to develop a
  //       flagging mechanism to tag materials for fission scaling.
  if (transient_options_.scale_fission_xs)
    for (const auto& xs : chi::trnsprt_xs_stack)
      if (not xs->IsFissionScaled())
        xs->ScaleFissionData(k_eff_);

  if (transient_options_.verbosity_level >= 1)
  {
    const double FR = ComputeFissionRate(false);
    char buff[200];
    snprintf(buff,200, " Initial Fission Rate FR=%12.6g", FR);
    chi::log.Log() << TextName() << buff;
  }

  //======================================== Compute auxiliary vectors
  fission_rate_local_.resize(grid_ptr_->local_cells.size(), 0.0);
  phi_prev_local_ = phi_old_local_;
  precursor_prev_local_ = precursor_new_local_;
  psi_prev_local_ = psi_new_local_;

  if (transient_options_.verbosity_level >= 0)
  {
    const double beta = ComputeBeta();
    char buff[200];
    snprintf(buff,200, " Beta=%.2f [pcm] reactivity=%.3f [$]",
            beta*1e5, (1.0- 1.0 / k_eff_) / beta);
    chi::log.Log() << TextName() << buff;
  }

  //================================================== Initialize source func
  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&DiscOrdTransientSolver::SetTransientSource, this, _1, _2, _3, _4);
}