#include "lbs_matrixaction_Ax.h"
#include "../Tools/ksp_data_context.h"
#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include "../../DiffusionSolver/Solver/diffusion_solver.h"

typedef chi_mesh::sweep_management::SweepScheduler MainSweepScheduler;
//###################################################################
/**Computes the action of the transport matrix on a vector.*/
int LinearBoltzmann::LBSMatrixAction_Ax(Mat matrix, Vec krylov_vector, Vec Ax)
{
  constexpr bool WITH_DELAYED_PSI = true;
  KSPDataContext* context;
  MatShellGetContext(matrix,&context);

  //Shorten some names
  LinearBoltzmann::Solver& solver = context->solver;
  SweepChunk& sweep_chunk = context->sweep_chunk;
  LBSGroupset& groupset  = context->groupset;
  MainSweepScheduler& sweepScheduler = context->sweepScheduler;
  SourceFlags& lhs_source_scope = context->lhs_scope;

  //============================================= Copy krylov vector into local
  solver.DisAssembleVector(groupset,
                           krylov_vector,
                           solver.phi_old_local.data(),WITH_DELAYED_PSI);

  //============================================= Setting the source using
  //                                             updated phi_old
  solver.SetSource(groupset, lhs_source_scope);

  //============================================= Sweeping the new source
  groupset.ZeroAngularFluxDataStructures();
  solver.phi_new_local.assign(solver.phi_new_local.size(),0.0);
  sweepScheduler.Sweep(sweep_chunk);

  //=================================================== Apply WGDSA
  if (groupset.apply_wgdsa)
  {
    solver.AssembleWGDSADeltaPhiVector(groupset,
                                       solver.phi_old_local.data(),
                                       solver.phi_new_local.data());
    ((chi_diffusion::Solver*)groupset.wgdsa_solver)->ExecuteS(true,false);
    solver.DisAssembleWGDSADeltaPhiVector(groupset,
                                          solver.phi_new_local.data());
  }
  if (groupset.apply_tgdsa)
  {
    solver.AssembleTGDSADeltaPhiVector(groupset,
                                       solver.phi_old_local.data(),
                                       solver.phi_new_local.data());
    ((chi_diffusion::Solver*)groupset.tgdsa_solver)->ExecuteS(true,false);
    solver.DisAssembleTGDSADeltaPhiVector(groupset,
                                          solver.phi_new_local.data());
  }


  solver.AssembleVector(groupset,
                        context->operating_vector,
                        solver.phi_new_local.data(),WITH_DELAYED_PSI);


  //============================================= Computing action
  VecWAXPY(Ax,-1.0,context->operating_vector,krylov_vector);

  return 0;
}