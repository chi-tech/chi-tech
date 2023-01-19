#include "lbts_transient_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Transient solver initialize routine.*/
void lbs::TransientSolver::Initialize()
{
  chi::log.Log() << "Initializing " << TextName() << ".";
  options.save_angular_flux = true;
  KEigenvalueSolver::Initialize();
  KEigenvalueSolver::Execute();

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
  if (transient_options.scale_fission_xs)
    for (const auto& xs : chi::trnsprt_xs_stack)
      if (!xs->is_fission_scaled)
        xs->ScaleFissionData(k_eff);

  if (transient_options.verbosity_level >= 1)
  {
    const double FR = ComputeFissionRate(false);
    char buff[200];
    sprintf(buff, " Initial Fission Rate FR=%12.6g", FR);
    chi::log.Log() << TextName() << buff;
  }

  //======================================== Compute auxiliary vectors
  fission_rate_local.resize(grid->local_cells.size(), 0.0);
  phi_prev_local = phi_old_local;
  precursor_prev_local = precursor_new_local;
  psi_prev_local = psi_new_local;

  if (transient_options.verbosity_level >= 0)
  {
    const double beta = ComputeBeta();
    char buff[200];
    sprintf(buff, " Beta=%.2f [pcm] reactivity=%.3f [$]",
            beta*1e5, (1.0-1.0/k_eff)/beta);
    chi::log.Log() << TextName() << buff;
  }

  //================================================== Initialize groupsets for
  //                                                   sweeping
  chi::log.Log() << "Initializing groupset sweeping data" << TextName() << ".";
  for (auto& groupset : groupsets)
  {
    ComputeSweepOrderings(groupset);
    InitFluxDataStructures(groupset);

    InitWGDSA(groupset);
    InitTGDSA(groupset);
  }

  chi::log.Log() << "Done initializing " << TextName() << ".";

  //================================================== Initialize source func
  using namespace std::placeholders;
  active_set_source_function =
    std::bind(&TransientSolver::SetTransientSource, this, _1, _2, _3);
}