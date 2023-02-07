#include "mip_wgs_context.h"

#include <petscksp.h>

#include "LBSMIPSteadyState/lbsmip_steady_solver.h"
#include "LBSMIPSteadyState/IterativeOperations/lbsmip_shell_operations.h"

#include "LBSSteadyState/Acceleration/diffusion_mip.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <iomanip>

#define sc_double static_cast<double>
#define PCShellPtr PetscErrorCode (*)(PC, Vec, Vec)

namespace lbs
{

template<>
void MIPWGSContext<Mat, Vec, KSP>::PreSetupCallback()
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
      << "********** Solving groupset " << groupset_.id
      << " with " << method_name << ".\n\n";
  }
}

template<>
void MIPWGSContext<Mat, Vec, KSP>::SetPreconditioner(KSP& solver)
{
  auto& ksp = solver;

  PC pc;
  KSPGetPC(ksp, &pc);

  if (groupset_.apply_wgdsa or groupset_.apply_tgdsa)
  {
    PCSetType(pc, PCSHELL);
    PCShellSetApply(pc, (PCShellPtr) MIP_TGDSA_PreConditionerMult);
    PCShellSetContext(pc, &(*this));
  }

  KSPSetPCSide(ksp,PC_LEFT);
  KSPSetUp(ksp);
}

template<>
std::pair<int64_t, int64_t> MIPWGSContext<Mat, Vec, KSP>::SystemSize()
{
  const size_t local_node_count = lbs_solver_.LocalNodeCount();
  const size_t globl_node_count = lbs_solver_.GlobalNodeCount();

  const size_t groupset_numgrps = groupset_.groups.size();
  const size_t local_size = local_node_count * groupset_numgrps;
  const size_t globl_size = globl_node_count * groupset_numgrps;

  return {static_cast<int64_t>(local_size),
          static_cast<int64_t>(globl_size)};
}

template<>
void MIPWGSContext<Mat, Vec, KSP>::ApplyInverseTransportOperator(int scope)
{
  auto& lbsmip_solver = dynamic_cast<MIPSteadyStateSolver&>(lbs_solver_);

  auto& mip_solver = *lbsmip_solver.gs_mip_solvers_[groupset_.id];

  std::vector<double> gs_q_moments_local_(SystemSize().first, 0.0);
  auto gs_phi_new_local_ = gs_q_moments_local_;

  lbsmip_solver.
    SetGSSTLvectorFromPrimarySTLvector(groupset_,
                                       gs_q_moments_local_,
                                       lbsmip_solver.QMomentsLocal());

  mip_solver.Assemble_b(gs_q_moments_local_);
  mip_solver.Solve(gs_phi_new_local_);

  lbsmip_solver.
    SetPrimarySTLvectorFromGSSTLvector(groupset_, gs_phi_new_local_,
                                       lbsmip_solver.PhiNewLocal());
}

template<>
void MIPWGSContext<Mat, Vec, KSP>::PostSolveCallback()
{
  lbs_solver_.GSScopedCopyPrimarySTLvectors(groupset_,
                                            lbs_solver_.PhiNewLocal(),
                                            lbs_solver_.PhiOldLocal());
}

}//namespace lbs