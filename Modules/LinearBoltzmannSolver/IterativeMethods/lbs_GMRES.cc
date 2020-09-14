#include "../lbs_linear_boltzmann_solver.h"

#include "ChiTimer/chi_timer.h"
#include "../Tools/kspmonitor_npt.h"
#include "../Tools/ksp_data_context.h"
#include "../IterativeOperations/lbs_matrixaction_Ax.h"

#include "../../DiffusionSolver/Solver/diffusion_solver.h"

#include "ChiPhysics/chi_physics.h"

#include <chi_log.h>
#include <chi_mpi.h>
extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;
extern ChiTimer chi_program_timer;

//###################################################################
/**Solves a groupset using GMRES.*/
void LinearBoltzmann::Solver::GMRES(int group_set_num, SweepChunk* sweep_chunk,
                                    MainSweepScheduler & sweepScheduler, bool log_info /* = true*/)
{
  constexpr bool WITH_DELAYED_PSI = true;
  if (log_info)
  {
    chi_log.Log(LOG_0)
      << "\n\n";
    chi_log.Log(LOG_0)
      << "********** Solving groupset " << group_set_num
      << " with GMRES.\n\n";
  }

  //================================================== Obtain groupset
  LBSGroupset* groupset = group_sets[group_set_num];
  int groupset_numgrps = groupset->groups.size();

  if (log_info)
    chi_log.Log(LOG_0)
      << "Quadrature number of angles: "
      << groupset->quadrature->abscissae.size() << "\n"
      << "Groups " << groupset->groups.front()->id << " "
      << groupset->groups.back()->id << "\n\n";

  //=================================================== Create Data context
  //                                                    available inside
  //                                                    Action
  KSPDataContext data_context;
  data_context.solver         = this;
  data_context.sweep_chunk    = sweep_chunk;
  data_context.group_set_num  = group_set_num;
  data_context.groupset       = groupset;
  data_context.sweepScheduler = &sweepScheduler;


  //=================================================== Create the matrix
  Mat A;
  auto num_ang_unknowns = groupset->angle_agg->GetNumberOfAngularUnknowns();
  int local_size = local_dof_count*num_moments*groupset_numgrps +
                   num_ang_unknowns.first;
  int globl_size = glob_dof_count*num_moments*groupset_numgrps +
                   num_ang_unknowns.second;
  MatCreateShell(PETSC_COMM_WORLD,local_size,
                                  local_size,
                                  globl_size,
                                  globl_size,
                                  &data_context,&A);
  groupset->angle_agg->ZeroIncomingDelayedPsi();

  //================================================== Set the action-operator
  MatShellSetOperation(A, MATOP_MULT, (void (*)(void)) LBSMatrixAction_Ax);

  //================================================== Initial vector assembly
  VecCreate(PETSC_COMM_WORLD,&phi_new);
  VecCreate(PETSC_COMM_WORLD,&phi_old);
  VecCreate(PETSC_COMM_WORLD,&q_fixed);

  VecSetSizes(phi_new,
              local_size,     //Local size
              globl_size);     //Global size
  VecSetType(phi_new,VECMPI);
  VecSet(phi_new,0.0);
  VecDuplicate(phi_new,&phi_old);
  VecDuplicate(phi_new,&q_fixed);
  VecDuplicate(phi_new,&data_context.x_temp);

  //================================================== Create Krylov Solver
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp,KSPGMRES);
  KSPSetOperators(ksp,A,A);
  data_context.krylov_solver = ksp;

  PC pc;
  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCNONE);

  KSPSetTolerances(ksp,1.e-50,
                   groupset->residual_tolerance,1.0e50,
                   groupset->max_iterations);
  KSPGMRESSetRestart(ksp,groupset->gmres_restart_intvl);
  KSPSetApplicationContext(ksp,&data_context);
  KSPSetConvergenceTest(ksp,&KSPConvergenceTestNPT,NULL,NULL);
  KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  KSPSetUp(ksp);

  //================================================== Compute b
  if (log_info)
    chi_log.Log(LOG_0) << chi_program_timer.GetTimeString() << " Computing b";

  SetSource(group_set_num,SourceFlags::USE_MATERIAL_SOURCE,
                          SourceFlags::SUPPRESS_PHI_OLD);
  sweep_chunk->SetDestinationPhi(&phi_new_local);

  phi_new_local.assign(phi_new_local.size(),0.0);
  sweepScheduler.Sweep(sweep_chunk);

  //=================================================== Apply DSA
  if (groupset->apply_wgdsa)
  {
    std::vector<double> phi_old_gmres(phi_old_local.size(),0.0);
    AssembleWGDSADeltaPhiVector(groupset, phi_old_gmres.data(), phi_new_local.data());
    ((chi_diffusion::Solver*)groupset->wgdsa_solver)->ExecuteS(true,false);
    DisAssembleWGDSADeltaPhiVector(groupset, phi_new_local.data());
  }
  if (groupset->apply_tgdsa)
  {
    std::vector<double> phi_old_gmres(phi_old_local.size(),0.0);
    AssembleTGDSADeltaPhiVector(groupset, phi_old_gmres.data(), phi_new_local.data());
    ((chi_diffusion::Solver*)groupset->tgdsa_solver)->ExecuteS(true,false);
    DisAssembleTGDSADeltaPhiVector(groupset, phi_new_local.data());
  }

  //=================================================== Assemble vectors
  AssembleVector(groupset,q_fixed,phi_new_local.data(),WITH_DELAYED_PSI);
  AssembleVector(groupset,phi_old,phi_old_local.data(),WITH_DELAYED_PSI);

  //=================================================== Retool for GMRES
  sweep_chunk->SetDestinationPhi(&phi_new_local);
  sweep_chunk->suppress_surface_src = true; //Action of Ax specific

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
    chi_log.Log(LOG_0)
      << chi_program_timer.GetTimeString() << " Starting iterations";
  KSPSolve(ksp,q_fixed,phi_new);
  //****************************************************

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  if (reason != KSP_CONVERGED_RTOL)
    chi_log.Log(LOG_0WARNING)
      << "GMRES solver failed. "
      << "Reason: " << chi_physics::GetPETScConvergedReasonstring(reason);


  DisAssembleVector(groupset, phi_new, phi_new_local.data(),WITH_DELAYED_PSI);
  DisAssembleVector(groupset, phi_new, phi_old_local.data(),WITH_DELAYED_PSI);

  //==================================================== Clean up
  KSPDestroy(&ksp);
  VecDestroy(&phi_new);
  VecDestroy(&phi_old);
  VecDestroy(&q_fixed);
  VecDestroy(&data_context.x_temp);
  MatDestroy(&A);



  double sweep_time = sweepScheduler.GetAverageSweepTime();
  double chunk_overhead_ratio = 1.0-sweepScheduler.GetAngleSetTimings()[2];
  double source_time=
    chi_log.ProcessEvent(source_event_tag,
                         ChiLog::EventOperation::AVERAGE_DURATION);
  size_t num_angles = groupset->quadrature->abscissae.size();
  long int num_unknowns = (long int)glob_dof_count*
                          (long int)num_angles*
                          (long int)groupset->groups.size();

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
      << sweep_time*1.0e9*chi_mpi.process_count/num_unknowns;
    chi_log.Log(LOG_0)
      << "        Number of unknowns per sweep:  " << num_unknowns;
    chi_log.Log(LOG_0)
      << "\n\n";
  }

  std::string sweep_log_file_name =
    std::string("GS_") + std::to_string(group_set_num) +
    std::string("_SweepLog_") + std::to_string(chi_mpi.location_id) +
    std::string(".log");
  groupset->PrintSweepInfoFile(sweepScheduler.sweep_event_tag,sweep_log_file_name);
}


