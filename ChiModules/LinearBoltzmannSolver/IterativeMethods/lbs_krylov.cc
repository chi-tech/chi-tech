#include "../lbs_linear_boltzmann_solver.h"

#include "../Tools/kspmonitor_npt.h"
#include "../Tools/ksp_data_context.h"
#include "../IterativeOperations/lbs_shell_operations.h"

#include "ChiMath/PETScUtils/petsc_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiTimer/chi_timer.h"
#include "LinearBoltzmannSolver/Groupset/lbs_groupset.h"

#define sc_double static_cast<double>
#define PCShellPtr PetscErrorCode (*)(PC, Vec, Vec)

//###################################################################
/**Solves a groupset using a general Krylov method.*/
bool lbs::SteadySolver::Krylov(LBSGroupset& groupset,
                               MainSweepScheduler& sweep_scheduler,
                               SourceFlags lhs_src_scope,
                               SourceFlags rhs_src_scope,
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
  const size_t local_size = local_node_count * num_moments * groupset_numgrps +
                            num_delayed_psi_info.first;
  const size_t globl_size = glob_node_count * num_moments * groupset_numgrps +
                            num_delayed_psi_info.second;
  const size_t num_angles = groupset.quadrature->abscissae.size();
  const size_t num_psi_global = glob_node_count *
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
  phi_new = chi_math::PETScUtils::CreateVector(static_cast<int64_t>(local_size),
                                               static_cast<int64_t>(globl_size));
  Vec x_temp;
  VecSet(phi_new,0.0);
  VecDuplicate(phi_new,&phi_old);
  VecDuplicate(phi_new,&q_fixed);
  VecDuplicate(phi_new,&x_temp);

  //=================================================== Create Data context
  //                                                    available inside
  //                                                    Action
  KSPDataContext data_context(*this, groupset, x_temp,
                              sweep_scheduler, lhs_src_scope);

  //=================================================== Create the matrix-shell
  Mat A;
  MatCreateShell(PETSC_COMM_WORLD,static_cast<int64_t>(local_size),
                                  static_cast<int64_t>(local_size),
                                  static_cast<int64_t>(globl_size),
                                  static_cast<int64_t>(globl_size),
                                  &data_context,&A);

  //================================================== Set the action-operator
  MatShellSetOperation(A, MATOP_MULT, (void (*)()) MatrixAction_Ax);

  //================================================== Create Krylov Solver
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
  auto init_q_moments_local = q_moments_local;
  SetSource(groupset, q_moments_local, rhs_src_scope);

  //Tool the sweep chunk
  auto& sweep_chunk = sweep_scheduler.GetSweepChunk();
  bool use_surface_source_flag = (rhs_src_scope & APPLY_MATERIAL_SOURCE) and
                                 (not options.use_src_moments);
  sweep_chunk.SetSurfaceSourceActiveFlag(use_surface_source_flag);
  sweep_chunk.ZeroIncomingDelayedPsi();

  //Sweep
  sweep_chunk.ZeroFluxDataStructures();
  sweep_scheduler.Sweep();

  //=================================================== Assemble vectors
  SetPETScVecFromSTLvector(groupset, q_fixed, phi_new_local, WITH_DELAYED_PSI);
  SetPETScVecFromSTLvector(groupset, phi_old, phi_old_local);

  PCApply(pc, q_fixed, x_temp);
  VecNorm(x_temp, NORM_2, &data_context.rhs_preconditioned_norm);

  //=================================================== Retool for GMRES
  sweep_chunk.SetSurfaceSourceActiveFlag(lhs_src_scope & APPLY_MATERIAL_SOURCE);
  sweep_chunk.SetDestinationPhi(phi_new_local);

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

  SetSTLvectorFromPETScVec(groupset, phi_new, phi_new_local, WITH_DELAYED_PSI);
  SetSTLvectorFromPETScVec(groupset, phi_new, phi_old_local, WITH_DELAYED_PSI);

  //================================================== Perform final sweep
  //                                                   with converged phi and
  //                                                   delayed psi dofs
  ZeroOutflowBalanceVars(groupset);

  sweep_chunk.SetDestinationPhi(phi_new_local);
  sweep_chunk.SetSurfaceSourceActiveFlag(rhs_src_scope & APPLY_MATERIAL_SOURCE);

  q_moments_local = init_q_moments_local;
  SetSource(groupset, q_moments_local, lhs_src_scope | rhs_src_scope);

  sweep_chunk.ZeroDestinationPhi();
  sweep_scheduler.Sweep();

  ScopedCopySTLvectors(groupset, phi_new_local, phi_old_local);

  //==================================================== Clean up
  KSPDestroy(&ksp);
  VecDestroy(&phi_new);
  VecDestroy(&phi_old);
  VecDestroy(&q_fixed);
  VecDestroy(&x_temp);
  MatDestroy(&A);

  //==================================================== Print solution info
  {
    double sweep_time = sweep_scheduler.GetAverageSweepTime();
    double chunk_overhead_ratio = 1.0 - sweep_scheduler.GetAngleSetTimings()[2];
    double source_time=
      chi::log.ProcessEvent(source_event_tag,
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


