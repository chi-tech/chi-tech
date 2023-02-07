#include "LBSMIPSteadyState/lbsmip_steady_solver.h"

#include "LBSSteadyState/Acceleration/diffusion_mip.h"

#include "LBSMIPSteadyState/Tools/lbsmip_kspmonitor_npt.h"
#include "LBSMIPSteadyState/Tools/lbsmip_ksp_data_context.h"
#include "LBSMIPSteadyState/IterativeOperations/lbsmip_shell_operations.h"
#include "LBSSteadyState/Groupset/lbs_groupset.h"

#include "ChiMath/PETScUtils/petsc_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiTimer/chi_timer.h"

#include <iomanip>

#define sc_double static_cast<double>
#define PCShellPtr PetscErrorCode (*)(PC, Vec, Vec)

//###################################################################
/**Solves a groupset using a general Krylov method.*/
bool lbs::MIPSteadyStateSolver::Krylov(LBSGroupset& groupset,
                                    SourceFlags lhs_src_scope,
                                    SourceFlags rhs_src_scope,
                                    const SetSourceFunction& set_source_function,
                                    bool log_info /* = true*/)
{
  constexpr bool NO_DELAYED_PSI = false;

  if (log_info)
  {
    std::string method_name;
    switch (groupset.iterative_method)
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
      << "\n\n";
    chi::log.Log()
      << "********** Solving groupset " << groupset.id
      << " with " << method_name << ".\n\n";
    chi::log.Log()
      << "Groups " << groupset.groups.front().id << " "
      << groupset.groups.back().id << "\n\n";
  }

  //================================================== Get groupset dof sizes
  const size_t groupset_numgrps = groupset.groups.size();
  const size_t local_size = local_node_count_ * num_moments_ * groupset_numgrps;
  const size_t globl_size = glob_node_count_ * num_moments_ * groupset_numgrps;

  //================================================== Create PETSc vectors
  auto phi_new = chi_math::PETScUtils::
  CreateVector(static_cast<int64_t>(local_size),
               static_cast<int64_t>(globl_size));

  Vec phi_old, q_fixed, x_temp;
  VecSet(phi_new,0.0);
  VecDuplicate(phi_new,&phi_old);
  VecDuplicate(phi_new,&q_fixed);
  VecDuplicate(phi_new,&x_temp);

  //=================================================== Create Data context
  //                                                    available inside
  //                                                    Action
  MIPKSPDataContext data_context(*this, groupset, x_temp,
                                 lhs_src_scope,
                                 set_source_function,
                                 phi_old_local_,
                                 q_moments_local_,
                                 phi_new_local_);

  //=================================================== Create the matrix-shell
  Mat A;
  MatCreateShell(PETSC_COMM_WORLD,static_cast<int64_t>(local_size),
                 static_cast<int64_t>(local_size),
                 static_cast<int64_t>(globl_size),
                 static_cast<int64_t>(globl_size),
                 &data_context,&A);

  //================================================== Set the action-operator
  MatShellSetOperation(A, MATOP_MULT, (void (*)()) MIPMatrixAction_Ax);

  //================================================== Create Krylov SteadyStateSolver
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  switch (groupset.iterative_method)
  {
    case IterativeMethod::KRYLOV_RICHARDSON: KSPSetType(ksp,KSPRICHARDSON); break;
    case IterativeMethod::KRYLOV_GMRES:      KSPSetType(ksp,KSPGMRES); break;
    case IterativeMethod::KRYLOV_BICGSTAB:   KSPSetType(ksp,KSPBCGS); break;
    default: KSPSetType(ksp,KSPGMRES); break;
  }
  KSPSetOperators(ksp,A,A);

  KSPSetTolerances(ksp,1.e-50,
                   groupset.residual_tolerance,1.0e50,
                   groupset.max_iterations);
  KSPSetApplicationContext(ksp, &data_context);
  KSPSetConvergenceTest(ksp, &MIPKSPConvergenceTest, nullptr, nullptr);
  KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);

  if (groupset.iterative_method == IterativeMethod::KRYLOV_GMRES)
  {
    KSPGMRESSetRestart(ksp, groupset.gmres_restart_intvl);
    KSPGMRESSetBreakdownTolerance(ksp, 1.0e6);
  }

  PC pc;
  KSPGetPC(ksp, &pc);

  if (groupset.apply_wgdsa or groupset.apply_tgdsa)
  {
    PCSetType(pc, PCSHELL);
    PCShellSetApply(pc, (PCShellPtr) MIP_TGDSA_PreConditionerMult);
    PCShellSetContext(pc, &data_context);
  }

  KSPSetPCSide(ksp,PC_LEFT);
  KSPSetUp(ksp);

  //================================================== Setup GS vectors
  std::vector<double> gs_q_moments_local_(local_size, 0.0);
  auto gs_phi_new_local_ = gs_q_moments_local_;

  //================================================== Compute b
  if (log_info)
  {
    chi::log.Log() << chi::program_timer.GetTimeString() << " Computing b";
  }

  //SetSource for RHS
  auto init_q_moments_local = q_moments_local_;
  set_source_function(groupset, q_moments_local_, PhiOldLocal(), rhs_src_scope);

  auto& mip_solver = *gs_mip_solvers_[groupset.id];

  SetGSSTLvectorFromPrimarySTLvector(groupset, gs_q_moments_local_, q_moments_local_);

  mip_solver.Assemble_b(gs_q_moments_local_);
  mip_solver.Solve(gs_phi_new_local_);

  SetPrimarySTLvectorFromGSSTLvector(groupset, gs_phi_new_local_, phi_new_local_);

  //=================================================== Assemble vectors
  SetGSPETScVecFromPrimarySTLvector(groupset, q_fixed, phi_new_local_, NO_DELAYED_PSI);
  SetGSPETScVecFromPrimarySTLvector(groupset, phi_old, phi_old_local_);

  PCApply(pc, q_fixed, x_temp);
  VecNorm(x_temp, NORM_2, &data_context.rhs_preconditioned_norm);

  double phi_old_norm=0.0;
  VecNorm(phi_old,NORM_2,&phi_old_norm);

  if (phi_old_norm > 1.0e-10)
  {
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    VecCopy(phi_old,phi_new);
    if (log_info)
      chi::log.Log() << "Using phi_old as initial guess.";
  }


  //**************** CALL GMRES SOLVE ******************
  if (log_info)
  {
    chi::log.Log()
      << chi::program_timer.GetTimeString() << " Starting iterations";
  }
  KSPSolve(ksp,q_fixed,phi_new);
  //****************************************************

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  if (reason != KSP_CONVERGED_RTOL)
    chi::log.Log0Warning()
      << "Krylov solver failed. "
      << "Reason: " << chi_physics::GetPETScConvergedReasonstring(reason);

  SetPrimarySTLvectorFromGSPETScVec(groupset, phi_new, phi_new_local_, NO_DELAYED_PSI);
  SetPrimarySTLvectorFromGSPETScVec(groupset, phi_new, phi_old_local_, NO_DELAYED_PSI);

  q_moments_local_ = init_q_moments_local;

  //==================================================== Clean up
  KSPDestroy(&ksp);
  VecDestroy(&phi_new);
  VecDestroy(&phi_old);
  VecDestroy(&q_fixed);
  VecDestroy(&x_temp);
  MatDestroy(&A);

  //==================================================== Print solution info
  {
    double source_time=
      chi::log.ProcessEvent(source_event_tag_,
                            chi_objects::ChiLog::EventOperation::AVERAGE_DURATION);

    if (log_info)
    {
      chi::log.Log()
        << "\n\n";
      chi::log.Log()
        << "        Set Src Time/sweep (s):        "
        << source_time;
      chi::log.Log()
        << "\n\n";
    }
  }

  return reason == KSP_CONVERGED_RTOL;
}

