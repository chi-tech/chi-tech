#include "../lbs_linear_boltzmann_solver.h"

#include "../Tools/kspmonitor_npt.h"
#include "../Tools/ksp_data_context.h"
#include "../IterativeOperations/lbs_matrixaction_Ax.h"

#include "DiffusionSolver/Solver/diffusion_solver.h"

#include "ChiPhysics/chi_physics.h"

#include "ChiMath/PETScUtils/petsc_utils.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Solves a groupset using GMRES.*/
bool LinearBoltzmann::Solver::GMRES(LBSGroupset& groupset,
                                    MainSweepScheduler& sweep_scheduler,
                                    SourceFlags lhs_src_scope,
                                    SourceFlags rhs_src_scope,
                                    bool log_info /* = true*/)
{
  constexpr bool WITH_DELAYED_PSI = true;
  if (log_info)
  {
    chi_log.Log(LOG_0)
      << "\n\n";
    chi_log.Log(LOG_0)
      << "********** Solving groupset " << groupset.id
      << " with GMRES.\n\n";
    chi_log.Log(LOG_0)
      << "Quadrature number of angles: "
      << groupset.quadrature->abscissae.size() << "\n"
      << "Groups " << groupset.groups.front().id << " "
      << groupset.groups.back().id << "\n\n";
  }

  //================================================== Get groupset dof sizes
  size_t groupset_numgrps = groupset.groups.size();
  auto num_delayed_ang_DOFs = groupset.angle_agg.GetNumDelayedAngularDOFs();
  size_t local_size = local_node_count * num_moments * groupset_numgrps +
                      num_delayed_ang_DOFs.first;
  size_t globl_size = glob_node_count * num_moments * groupset_numgrps +
                      num_delayed_ang_DOFs.second;

  if (log_info)
    chi_log.Log(LOG_0)
      << "Number of angular unknowns: " << num_delayed_ang_DOFs.second;

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
  MatShellSetOperation(A, MATOP_MULT, (void (*)()) LBSMatrixAction_Ax);

  //================================================== Create Krylov Solver
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp,KSPGMRES);
  KSPSetOperators(ksp,A,A);

  KSPSetTolerances(ksp,1.e-50,
                   groupset.residual_tolerance,1.0e50,
                   groupset.max_iterations);
  KSPGMRESSetRestart(ksp, groupset.gmres_restart_intvl);
  KSPSetApplicationContext(ksp, &data_context);
  KSPSetConvergenceTest(ksp, &KSPConvergenceTestNPT, nullptr, nullptr);
  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  KSPSetUp(ksp);

  //================================================== Compute b
  if (log_info)
  {
    chi_log.Log(LOG_0) << chi_program_timer.GetTimeString() << " Computing b";
  }

  //SetSource for RHS
  auto init_q_moments_local = q_moments_local;
  SetSource(groupset, q_moments_local, rhs_src_scope);

  //Tool the sweep chunk
  auto& sweep_chunk = sweep_scheduler.sweep_chunk;
  bool use_surface_source_flag = (rhs_src_scope & APPLY_MATERIAL_SOURCE) and
                                 (not options.use_src_moments);
  sweep_chunk.SetSurfaceSourceActiveFlag(use_surface_source_flag);
  sweep_chunk.SetDestinationPhi(phi_new_local);
  sweep_chunk.ZeroIncomingDelayedPsi();

  //Sweep
  sweep_chunk.ZeroFluxDataStructures();
  sweep_scheduler.Sweep();

  //=================================================== Apply DSA
  if (groupset.apply_wgdsa)
  {
    AssembleWGDSADeltaPhiVector(groupset, phi_old_local.data(), phi_new_local.data());
    ((chi_diffusion::Solver*)groupset.wgdsa_solver)->ExecuteS(true,false);
    DisAssembleWGDSADeltaPhiVector(groupset, phi_new_local.data());
  }
  if (groupset.apply_tgdsa)
  {
    AssembleTGDSADeltaPhiVector(groupset, phi_old_local.data(), phi_new_local.data());
    ((chi_diffusion::Solver*)groupset.tgdsa_solver)->ExecuteS(true,false);
    DisAssembleTGDSADeltaPhiVector(groupset, phi_new_local.data());
  }

  //=================================================== Assemble vectors
  SetPETScVecFromSTLvector(groupset, q_fixed, phi_new_local, WITH_DELAYED_PSI);
  SetPETScVecFromSTLvector(groupset, phi_old, phi_old_local, WITH_DELAYED_PSI);

  //=================================================== Retool for GMRES
  sweep_chunk.SetSurfaceSourceActiveFlag(lhs_src_scope & APPLY_MATERIAL_SOURCE);
  sweep_chunk.SetDestinationPhi(phi_new_local);

  double phi_old_norm=0.0;
  VecNorm(phi_old,NORM_2,&phi_old_norm);

  if (phi_old_norm > 1.0e-10)
  {
    VecCopy(phi_old,phi_new);
    if (log_info)
      chi_log.Log(LOG_0) << "Using phi_old as initial guess.";
  }


  //**************** CALL GMRES SOLVE ******************
  if (log_info)
  {
    chi_log.Log(LOG_0)
      << chi_program_timer.GetTimeString() << " Starting iterations";
  }
  KSPSolve(ksp,q_fixed,phi_new);
  //****************************************************

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  if (reason != KSP_CONVERGED_RTOL)
    chi_log.Log(LOG_0WARNING)
      << "GMRES solver failed. "
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
      chi_log.ProcessEvent(source_event_tag,
                           ChiLog::EventOperation::AVERAGE_DURATION);
    size_t num_angles = groupset.quadrature->abscissae.size();
    size_t num_unknowns = glob_node_count *
                          num_angles *
                          groupset.groups.size();

    if (log_info)
    {
      chi_log.Log(LOG_0)
        << "\n\n";
      chi_log.Log(LOG_0)
        << "        Set Src Time/sweep (s):        "
        << source_time;
      chi_log.Log(LOG_0)
        << "        Average sweep time (s):        "
        << sweep_time;
      chi_log.Log(LOG_0)
        << "        Chunk-Overhead-Ratio  :        "
        << chunk_overhead_ratio;
      chi_log.Log(LOG_0)
        << "        Sweep Time/Unknown (ns):       "
        << sweep_time*1.0e9*chi_mpi.process_count/
           static_cast<double>(num_unknowns);
      chi_log.Log(LOG_0)
        << "        Number of unknowns per sweep:  " << num_unknowns;
      chi_log.Log(LOG_0)
        << "\n\n";

      std::string sweep_log_file_name =
          std::string("GS_") + std::to_string(groupset.id) +
          std::string("_SweepLog_") + std::to_string(chi_mpi.location_id) +
          std::string(".log");
      groupset.PrintSweepInfoFile(sweep_scheduler.sweep_event_tag,
                                  sweep_log_file_name);
    }
  }

  return reason == KSP_CONVERGED_RTOL;
}


