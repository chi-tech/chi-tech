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

  for (double& phi_value : phi_old_local)
    phi_value /= k_eff;

  ComputePrecursors();

  if (transient_options.verbosity_level >= 1)
  {
    const double FR = ComputeFissionRate(false);
    char buff[200];
    snprintf(buff,200, " Initial Fission Rate FR=%12.6g", FR);
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
    snprintf(buff,200, " Beta=%.2f [pcm] reactivity=%.3f [$]",
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