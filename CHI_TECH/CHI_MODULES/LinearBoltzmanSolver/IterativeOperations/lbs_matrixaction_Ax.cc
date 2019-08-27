#include "lbs_matrixaction_Ax.h"
#include "../Tools/ksp_data_context.h"
#include "CHI_MESH/CHI_SWEEP/chi_sweepscheduler.h"

#include <CHI_MODULES/CHI_DIFFUSION/Solver/diffusion_solver.h>

typedef chi_mesh::SweepManagement::SweepScheduler MainSweepScheduler;
//###################################################################
/**Computes the action of the transport matrix on a vector.*/
int NPTMatrixAction_Ax(Mat matrix, Vec krylov_vector, Vec Ax)
{
  KSP_DATA_CONTEXT* context;
  MatShellGetContext(matrix,&context);

  LinearBoltzmanSolver* solver = context->solver;
  SweepChunk* sweep_chunk = context->sweep_chunk;
  LBS_GROUPSET* groupset  = context->groupset;
  MainSweepScheduler* sweepScheduler = context->sweepScheduler;

  //============================================= Copy krylov vector into local
  solver->DisAssembleVector(groupset,
                            krylov_vector,
                            solver->phi_old_local.data());

  //============================================= Setting the source using
  //                                             updated phi_old
  solver->SetSource(context->group_set_num,USE_DLINV_SOURCE);

  //============================================= Sweeping the new source
  solver->phi_new_local.assign(solver->phi_new_local.size(),0.0);
  sweepScheduler->Sweep(sweep_chunk);

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