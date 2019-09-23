#include "lbs_matrixaction_Ax.h"
#include "../Tools/ksp_data_context.h"
#include "ChiMesh/SweepUtilities/chi_sweepscheduler.h"

#include "../../DiffusionSolver/Solver/diffusion_solver.h"

#include <ChiTimer/chi_timer.h>

extern ChiTimer chi_program_timer;

typedef chi_mesh::SweepManagement::SweepScheduler MainSweepScheduler;
//###################################################################
/**Computes the action of the transport matrix on a vector.*/
int NPTMatrixAction_Ax_Cycles(Mat matrix, Vec krylov_vector, Vec Ax)
{
  KSPDataContext* context;
  MatShellGetContext(matrix,&context);

  LinearBoltzmanSolver* solver = context->solver;
  SweepChunk* sweep_chunk = context->sweep_chunk;
  LBSGroupset* groupset  = context->groupset;
  MainSweepScheduler* sweepScheduler = context->sweepScheduler;

  //============================================= Copy krylov vector into local
  solver->DisAssembleVector(groupset,
                            krylov_vector,
                            solver->phi_old_local.data());

  //============================================= Setting the source using
  //                                             updated phi_old
  solver->SetSource(context->group_set_num,USE_DLINV_SOURCE);

  //============================================= Sweeping the new source
//  solver->phi_new_local.assign(solver->phi_new_local.size(),0.0);
//  sweepScheduler->Sweep(sweep_chunk);
  groupset->angle_agg->ResetDelayedPsi();
  double new_norm = 0.0;
  double prev_norm = 1.0;
  bool cycles_converged = false;
  for (int k=0; k<15; k++)
  {
    solver->phi_new_local.assign(solver->phi_new_local.size(),0.0); //Ensure phi_new=0.0
    sweepScheduler->Sweep(sweep_chunk);
    new_norm = groupset->angle_agg->GetDelayedPsiNorm();

    double rel_change = std::fabs(1.0-new_norm/prev_norm);
    prev_norm = new_norm;
//    if (rel_change<std::max(1.0e-2,1.0e-10))
//      cycles_converged = true;

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

    if (cycles_converged)
      break;
  }



  //=================================================== Apply WGDSA
  if (groupset->apply_wgdsa)
  {
    solver->AssembleWGDSADeltaPhiVector(groupset,
                                        solver->phi_old_local.data(),
                                        solver->phi_new_local.data());
    ((chi_diffusion::Solver*)groupset->wgdsa_solver)->ExecuteS(true,false);
    solver->DisAssembleWGDSADeltaPhiVector(groupset,
                                           solver->phi_new_local.data());
  }
  if (groupset->apply_tgdsa)
  {
    solver->AssembleTGDSADeltaPhiVector(groupset,
                                        solver->phi_old_local.data(),
                                        solver->phi_new_local.data());
    ((chi_diffusion::Solver*)groupset->tgdsa_solver)->ExecuteS(true,false);
    solver->DisAssembleTGDSADeltaPhiVector(groupset,
                                           solver->phi_new_local.data());
  }


  solver->AssembleVector(groupset,
                         context->x_temp,
                         solver->phi_new_local.data());


  //============================================= Computing action
  VecWAXPY(Ax,-1.0,context->x_temp,krylov_vector);

  return 0;
}