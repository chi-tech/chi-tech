#include "../lbs_linear_boltzman_solver.h"

#include "ChiMesh/SweepUtilities/chi_sweepscheduler.h"
#include "../SweepChunks/lbs_sweepchunk_pwl_polyhedron.h"

#include "ChiTimer/chi_timer.h"
#include "../Tools/kspmonitor_npt.h"
#include "../Tools/ksp_data_context.h"
#include "../IterativeOperations/lbs_matrixaction_Ax.h"

#include "../../DiffusionSolver/Solver/diffusion_solver.h"

#include <chi_log.h>

extern ChiLog chi_log;

extern double chi_global_timings[20];
extern ChiTimer chi_program_timer;

typedef chi_mesh::SweepManagement::SweepChunk SweepChunk;
typedef chi_mesh::SweepManagement::SweepScheduler MainSweepScheduler;

//###################################################################
/**Solves a groupset using GMRES.*/
void LinearBoltzman::Solver::GMRES_Cycles(int group_set_num)
{
  chi_log.Log(LOG_0)
    << "\n\n";
  chi_log.Log(LOG_0)
    << "********** Solving groupset " << group_set_num
    << " with GMRES.\n\n";

  //================================================== Obtain groupset
  LBSGroupset* groupset = group_sets[group_set_num];
  int groupset_numgrps = groupset->groups.size();
  chi_log.Log(LOG_0)
    << "Quadrature number of angles: "
    << groupset->quadrature->abscissae.size() << "\n"
    << "Number of azimuthal angles : "
    << groupset->quadrature->azimu_ang.size() << "\n"
    << "Number of polar angles     : "
    << groupset->quadrature->polar_ang.size() << "\n\n";

  //================================================== Setting up required
  //                                                   sweep chunks
  SweepChunk* sweep_chunk = SetSweepChunk(group_set_num);
  MainSweepScheduler sweepScheduler(DEPTH_OF_GRAPH,
                                    groupset->angle_agg);

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
  int local_size = local_dof_count*num_moments*groupset_numgrps;
  int globl_size = glob_dof_count*num_moments*groupset_numgrps;
  MatCreateShell(PETSC_COMM_WORLD,local_size,
                                  local_size,
                                  globl_size,
                                  globl_size,
                                  &data_context,&A);

  //================================================== Set the action-operator
  MatShellSetOperation(A, MATOP_MULT,(void (*)(void)) NPTMatrixAction_Ax_Cycles);

  //================================================== Initial vector assembly
  phi_new_local.assign(phi_new_local.size(),0.0);
  phi_oldcycle_local.resize(phi_new_local.size(),0.0);

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
  chi_log.Log(LOG_0) << "Computing b";
  SetSource(group_set_num,SourceFlags::USE_MATERIAL_SOURCE,SourceFlags::SUPPRESS_PHI_OLD);
  sweep_chunk->SetDestinationPhi(&phi_new_local);

  double new_norm = 0.0;
  double prev_norm = 1.0;
  bool cycles_converged = false;
  bool sweep_once_more = false;
  for (int k=0; k<20; k++)
  {
    phi_new_local.assign(phi_new_local.size(),0.0); //Ensure phi_new=0.0
    sweepScheduler.Sweep(sweep_chunk);
    new_norm = groupset->angle_agg->GetDelayedPsiNorm();

    double rel_change = 0.0;
    if (prev_norm > 1.0e-10)
      rel_change = std::fabs(1.0-new_norm/prev_norm);
    prev_norm = new_norm;
    if (new_norm<std::max(1.0e-8,1.0e-10))
      cycles_converged = true;

    std::string offset;
    if (groupset->apply_wgdsa || groupset->apply_tgdsa)
      offset = std::string("    ");

    std::stringstream iter_info;
    iter_info
        << chi_program_timer.GetTimeString() << " "
        << offset
        << "Cyclic iteration " << std::setw(5) << k
        << " Point-wise change " << std::setw(14) << rel_change << " " << new_norm;

    if (cycles_converged)
      iter_info << " CONVERGED\n";

    chi_log.Log(LOG_0) << iter_info.str();

    if (cycles_converged and (!sweep_once_more))
      sweep_once_more = true;
    else if (cycles_converged)
      break;
  }



  //=================================================== Apply WGDSA
  if (groupset->apply_wgdsa)
  {
    AssembleWGDSADeltaPhiVector(groupset, phi_old_local.data(), phi_new_local.data());
    ((chi_diffusion::Solver*)groupset->wgdsa_solver)->ExecuteS(true,false);
    DisAssembleWGDSADeltaPhiVector(groupset, phi_new_local.data());
  }
  if (groupset->apply_tgdsa)
  {
    AssembleTGDSADeltaPhiVector(groupset, phi_old_local.data(), phi_new_local.data());
    ((chi_diffusion::Solver*)groupset->tgdsa_solver)->ExecuteS(true,false);
    DisAssembleTGDSADeltaPhiVector(groupset, phi_new_local.data());
  }

  //=================================================== Assemble vectors
  AssembleVector(groupset,q_fixed,phi_new_local.data());
  AssembleVector(groupset,phi_old,phi_old_local.data());

  //=================================================== Retool for GMRES
  sweep_chunk->SetDestinationPhi(&phi_new_local);
  sweep_chunk->suppress_surface_src = true; //Action of Ax specific

  VecCopy(phi_old,phi_new);

  //**************** CALL GMRES SOLVE ******************
  chi_log.Log(LOG_0) << "Starting iterations";
  KSPSolve(ksp,q_fixed,phi_new);
//  AssembleVector(groupset,phi_new,phi_new_local.data());
  //****************************************************

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp,&reason);
  if (reason<0)
    chi_log.Log(LOG_0WARNING) << "GMRES solver failed.";


  DisAssembleVector(groupset, phi_new, phi_new_local.data());
  DisAssembleVector(groupset, phi_new, phi_old_local.data());

  //==================================================== Clean up
  KSPDestroy(&ksp);
  VecDestroy(&phi_new);
  VecDestroy(&phi_old);
  VecDestroy(&q_fixed);



  double sweep_time = chi_global_timings[16]/chi_global_timings[17];
  double source_time= chi_global_timings[18]/chi_global_timings[19];
  size_t num_angles = groupset->quadrature->abscissae.size();
  long int num_unknowns = (long int)glob_dof_count*
                          (long int)num_angles*
                          (long int)groupset->groups.size();
  chi_log.Log(LOG_0)
    << "\n\n";
  chi_log.Log(LOG_0)
    << "        Set Src Time/sweep (s):    "
    << source_time*1.0e-3;
  chi_log.Log(LOG_0)
    << "        Sweep Time/Unknown (ns):       "
    << sweep_time*1.0e9*chi_mpi.process_count/num_unknowns;
  chi_log.Log(LOG_0)
    << "        Number of unknowns per sweep:  " << num_unknowns;
  chi_log.Log(LOG_0)
    << "\n\n";
}


