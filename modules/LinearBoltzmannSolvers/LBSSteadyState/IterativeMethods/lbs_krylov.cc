#include "LBSSteadyState/lbs_linear_boltzmann_solver.h"

#include "LBSSteadyState/Tools/kspmonitor_npt.h"
#include "LBSSteadyState/Tools/ksp_data_context.h"
#include "LBSSteadyState/IterativeOperations/lbs_shell_operations.h"
#include "LBSSteadyState/Groupset/lbs_groupset.h"

#include "ChiMath/PETScUtils/petsc_utils.h"

#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiTimer/chi_timer.h"

#include <iomanip>

#define sc_double static_cast<double>
#define PCShellPtr PetscErrorCode (*)(PC, Vec, Vec)

//###################################################################
/**Solves a groupset using a general Krylov method.*/
bool lbs::SteadyStateSolver::
  Krylov(LBSGroupset& groupset,
         chi_mesh::sweep_management::SweepScheduler& sweep_scheduler,
         SourceFlags lhs_src_scope,
         SourceFlags rhs_src_scope,
         const SetSourceFunction& set_source_function,
         bool log_info /* = true*/)
{
  constexpr bool WITH_DELAYED_PSI = true;

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
      << "Quadrature number of angles: "
      << groupset.quadrature->abscissae.size() << "\n"
      << "Groups " << groupset.groups.front().id << " "
      << groupset.groups.back().id << "\n\n";
  }

  //================================================== Get groupset dof sizes
  const size_t groupset_numgrps = groupset.groups.size();
  const auto num_delayed_psi_info = groupset.angle_agg.GetNumDelayedAngularDOFs();
  const size_t local_size = local_node_count_ * num_moments_ * groupset_numgrps +
                            num_delayed_psi_info.first;
  const size_t globl_size = glob_node_count_ * num_moments_ * groupset_numgrps +
                            num_delayed_psi_info.second;
  const size_t num_angles = groupset.quadrature->abscissae.size();
  const size_t num_psi_global = glob_node_count_ *
                                num_angles *
                                groupset.groups.size();
  const size_t num_delayed_psi_globl = num_delayed_psi_info.second;

  if (log_info)
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


  //================================================== Create PETSc vectors
  auto phi_new = chi_math::PETScUtils::
    CreateVector(static_cast<int64_t>(local_size),
                 static_cast<int64_t>(globl_size));

  Vec q_fixed;
  VecSet(phi_new,0.0);
  VecDuplicate(phi_new,&q_fixed);

  //=================================================== Create Data context
  //                                                    available inside
  //                                                    Action
  KSPDataContext data_context(*this, groupset,
                              sweep_scheduler, lhs_src_scope,
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
  MatShellSetOperation(A, MATOP_MULT, (void (*)()) MatrixAction_Ax);

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
  KSPSetConvergenceTest(ksp, &KSPConvergenceTest, nullptr, nullptr);
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
    PCShellSetApply(pc, (PCShellPtr) WGDSA_TGDSA_PreConditionerMult);
    PCShellSetContext(pc, &data_context);
  }

  KSPSetPCSide(ksp,PC_LEFT);
  KSPSetUp(ksp);

  //================================================== Compute b
  if (log_info)
  {
    chi::log.Log() << chi::program_timer.GetTimeString() << " Computing b";
  }

  //SetSource for RHS
  auto init_q_moments_local = q_moments_local_;
  set_source_function(groupset, q_moments_local_, PhiOldLocal(), rhs_src_scope);

  //Tool the sweep chunk
  auto& sweep_chunk = sweep_scheduler.GetSweepChunk();
  bool use_surface_source_flag = (rhs_src_scope & APPLY_FIXED_SOURCES) and
                                 (not options_.use_src_moments);
  sweep_chunk.SetSurfaceSourceActiveFlag(use_surface_source_flag);
  sweep_chunk.ZeroIncomingDelayedPsi();

  //Sweep
  sweep_chunk.ZeroFluxDataStructures();
  sweep_scheduler.Sweep();

  //=================================================== Assemble vectors
  SetGSPETScVecFromPrimarySTLvector(groupset, q_fixed, phi_new_local_, WITH_DELAYED_PSI);
  SetGSPETScVecFromPrimarySTLvector(groupset, phi_new, phi_old_local_);

  Vec x_temp;
  VecDuplicate(q_fixed, &x_temp);
  PCApply(pc, q_fixed, x_temp);
  VecNorm(x_temp, NORM_2, &data_context.rhs_preconditioned_norm);
  VecDestroy(&x_temp);

  //=================================================== Retool for GMRES
  sweep_chunk.SetSurfaceSourceActiveFlag(lhs_src_scope & APPLY_FIXED_SOURCES);
  sweep_chunk.SetDestinationPhi(phi_new_local_);

  double phi_old_norm=0.0;
  VecNorm(phi_new,NORM_2,&phi_old_norm);

  if (phi_old_norm > 1.0e-10)
  {
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
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

  SetPrimarySTLvectorFromGSPETScVec(groupset, phi_new, phi_new_local_, WITH_DELAYED_PSI);
  SetPrimarySTLvectorFromGSPETScVec(groupset, phi_new, phi_old_local_, WITH_DELAYED_PSI);

  //================================================== Perform final sweep
  //                                                   with converged phi and
  //                                                   delayed psi dofs
  ZeroOutflowBalanceVars(groupset);

  sweep_chunk.SetDestinationPhi(phi_new_local_);
  sweep_chunk.SetSurfaceSourceActiveFlag(rhs_src_scope & APPLY_FIXED_SOURCES);

  q_moments_local_ = init_q_moments_local;
  set_source_function(groupset, q_moments_local_,
                      PhiOldLocal(),lhs_src_scope | rhs_src_scope);

  sweep_chunk.ZeroDestinationPhi();
  sweep_scheduler.Sweep();

  GSScopedCopyPrimarySTLvectors(groupset, phi_new_local_, phi_old_local_);

  //==================================================== Clean up
  KSPDestroy(&ksp);
  VecDestroy(&phi_new);
  VecDestroy(&q_fixed);
  VecDestroy(&x_temp);
  MatDestroy(&A);

  //==================================================== Print solution info
  {
    double sweep_time = sweep_scheduler.GetAverageSweepTime();
    double chunk_overhead_ratio = 1.0 - sweep_scheduler.GetAngleSetTimings()[2];
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
        << "        Average sweep time (s):        "
        << sweep_time;
      chi::log.Log()
        << "        Chunk-Overhead-Ratio  :        "
        << chunk_overhead_ratio;
      chi::log.Log()
        << "        Sweep Time/Unknown (ns):       "
        << sweep_time*1.0e9*chi::mpi.process_count/
           sc_double(num_psi_global);
      chi::log.Log()
        << "        Number of unknowns per sweep:  " << num_psi_global;
      chi::log.Log()
        << "\n\n";

      std::string sweep_log_file_name =
          std::string("GS_") + std::to_string(groupset.id) +
          std::string("_SweepLog_") + std::to_string(chi::mpi.location_id) +
          std::string(".log");
      groupset.PrintSweepInfoFile(sweep_scheduler.sweep_event_tag,
                                  sweep_log_file_name);
    }
  }

  return reason == KSP_CONVERGED_RTOL;
}


