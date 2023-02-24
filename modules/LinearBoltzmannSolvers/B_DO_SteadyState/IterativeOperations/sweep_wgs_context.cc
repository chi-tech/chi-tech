#include "sweep_wgs_context.h"

#include <petscksp.h>

#include "LinearBoltzmannSolvers/B_DO_SteadyState/lbs_DO_steady_state.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Preconditioning/lbs_shell_operations.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <iomanip>

#define sc_double static_cast<double>
#define PCShellPtr PetscErrorCode (*)(PC, Vec, Vec)

namespace lbs
{
/**General print out of information.*/
template<>
void SweepWGSContext<Mat, Vec, KSP>::PreSetupCallback()
{
  if (log_info_)
  {
    std::string method_name;
    switch (groupset_.iterative_method)
    {
      case IterativeMethod::KRYLOV_RICHARDSON:
        method_name = "KRYLOV_RICHARDSON"; break;
      case IterativeMethod::KRYLOV_GMRES:
        method_name = "KRYLOV_GMRES"; break;
      case IterativeMethod::KRYLOV_BICGSTAB:
        method_name = "KRYLOV_BICGSTAB"; break;
      default: method_name = "KRYLOV_GMRES";
    }
    chi::log.Log()
      << "\n\n"
      << "********** Solving groupset " << groupset_.id_
      << " with " << method_name << ".\n\n"
      << "Quadrature number of angles: "
      << groupset_.quadrature->abscissae.size() << "\n"
      << "Groups " << groupset_.groups_.front().id_ << " "
      << groupset_.groups_.back().id_ << "\n\n";
  }
}

/**Sets the preconditioner application function.*/
template<>
void SweepWGSContext<Mat, Vec, KSP>::SetPreconditioner(KSP& solver)
{
  auto& ksp = solver;

  PC pc;
  KSPGetPC(ksp, &pc);

  if (groupset_.apply_wgdsa or groupset_.apply_tgdsa)
  {
    PCSetType(pc, PCSHELL);
    PCShellSetApply(pc, (PCShellPtr) WGDSA_TGDSA_PreConditionerMult);
    PCShellSetContext(pc, &(*this));
  }

  KSPSetPCSide(ksp,PC_LEFT);
  KSPSetUp(ksp);
}

/**For sweeping we add lagged angular fluxes to the size of the vectors.*/
template<>
std::pair<int64_t, int64_t> SweepWGSContext<Mat, Vec, KSP>::SystemSize()
{
  const size_t local_node_count = lbs_solver_.LocalNodeCount();
  const size_t globl_node_count = lbs_solver_.GlobalNodeCount();
  const size_t num_moments      = lbs_solver_.NumMoments();

  const size_t groupset_numgrps = groupset_.groups_.size();
  const auto num_delayed_psi_info = groupset_.angle_agg.GetNumDelayedAngularDOFs();
  const size_t local_size = local_node_count * num_moments * groupset_numgrps +
                            num_delayed_psi_info.first;
  const size_t globl_size = globl_node_count * num_moments * groupset_numgrps +
                            num_delayed_psi_info.second;
  const size_t num_angles = groupset_.quadrature->abscissae.size();
  const size_t num_psi_global = globl_node_count *
                                num_angles *
                                groupset_.groups_.size();
  const size_t num_delayed_psi_globl = num_delayed_psi_info.second;

  if (log_info_)
  {
    chi::log.Log()
      << "Total number of angular unknowns: "
      << num_psi_global
      << "\n"
      << "Number of lagged angular unknowns: "
      << num_delayed_psi_globl << "("
      << std::setprecision(2)
      << sc_double(num_delayed_psi_globl)*100 / sc_double(num_psi_global)
      << "%)";
  }

  return {static_cast<int64_t>(local_size),
          static_cast<int64_t>(globl_size)};
}

/**With a right-hand side built. This routine applies the inverse
 * of the transport operator to this right-hand side.*/
template<>
void SweepWGSContext<Mat, Vec, KSP>::ApplyInverseTransportOperator(int scope)
{
  const bool use_surface_source_flag =
    (scope & APPLY_FIXED_SOURCES) and
    (not lbs_solver_.Options().use_src_moments);

  auto& sweep_chunk = sweep_scheduler_.GetSweepChunk();

  sweep_chunk.SetSurfaceSourceActiveFlag(use_surface_source_flag);

  if (scope & APPLY_FIXED_SOURCES)
    sweep_chunk.ZeroIncomingDelayedPsi();

  //Sweep
  sweep_chunk.ZeroFluxDataStructures();
  sweep_scheduler_.Sweep();
}

/**This method implements an additional sweep for two reasons:
 * The first is to compute balance parameters, and the second
 * is to allow for the calculation of proper angular fluxes. The
 * latter is needed because some krylov methods do not necessarily
 * provide the true angular flux at each iteration.*/
template<>
void SweepWGSContext<Mat, Vec, KSP>::PostSolveCallback()
{
  auto& sweep_chunk = sweep_scheduler_.GetSweepChunk();

  //================================================== Perform final sweep
  //                                                   with converged phi and
  //                                                   delayed psi dofs
  lbs_ss_solver_.ZeroOutflowBalanceVars(groupset_);

  sweep_chunk.SetDestinationPhi(lbs_solver_.PhiNewLocal());
  sweep_chunk.SetSurfaceSourceActiveFlag(rhs_src_scope_ & APPLY_FIXED_SOURCES);

  set_source_function_(groupset_, lbs_solver_.QMomentsLocal(),
                       lbs_solver_.PhiOldLocal(),
                       lhs_src_scope_ | rhs_src_scope_);

  sweep_chunk.ZeroDestinationPhi();
  sweep_scheduler_.Sweep();

  lbs_solver_.GSScopedCopyPrimarySTLvectors(groupset_,
                                            PhiSTLOption::PHI_NEW,
                                            PhiSTLOption::PHI_OLD);
}

}//namespace lbs