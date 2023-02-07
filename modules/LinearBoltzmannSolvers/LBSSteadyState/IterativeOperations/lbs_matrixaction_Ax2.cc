#include "lbs_shell_operations.h"
#include "LBSSteadyState/Tools/ksp_data_context.h"
#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include "LBSSteadyState/Groupset/lbs_groupset.h"

typedef chi_mesh::sweep_management::SweepScheduler MainSweepScheduler;
//###################################################################
/**Computes the action of the transport matrix on a vector.*/
int lbs::MatrixAction_Ax(Mat matrix, Vec krylov_vector, Vec Av)
{
  constexpr bool WITH_DELAYED_PSI = true;
  KSPDataContext* context;
  MatShellGetContext(matrix,&context);

  //Shorten some names
  lbs::SteadyStateSolver& solver = context->solver;
  LBSGroupset& groupset  = context->groupset;
  MainSweepScheduler& sweep_scheduler = context->sweep_scheduler;
  auto& sweep_chunk = context->sweep_scheduler.GetSweepChunk();
  SourceFlags& lhs_source_scope = context->lhs_scope;
  auto& set_source_function = context->set_source_function;

  //============================================= Copy krylov vector into local
  solver.SetPrimarySTLvectorFromGSPETScVec(groupset,
                                           krylov_vector,
                                           context->phi_old_local, WITH_DELAYED_PSI);

  //============================================= Setting the source using
  //                                              updated phi_old
  context->q_moments_local.assign(context->q_moments_local.size(), 0.0);
  set_source_function(groupset,
                      context->q_moments_local,
                      solver.PhiOldLocal(),
                      lhs_source_scope);

  //============================================= Sweeping the new source
  sweep_chunk.ZeroFluxDataStructures();
  sweep_scheduler.Sweep();

  //=================================================== Copy local into
  //                                                    operating vector
  // We copy the STL data to the operating vector
  // petsc_phi_delta first because its already sized.
  // pc_output is not necessarily initialized yet.
  solver.SetGSPETScVecFromPrimarySTLvector(groupset,
                                           Av,
                                           context->phi_new_local, WITH_DELAYED_PSI);

  //============================================= Computing action
  // A  = [I - DLinvMS]
  // Av = [I - DLinvMS]v
  //    = v - DLinvMSv
  VecAYPX(Av, -1.0, krylov_vector);

  return 0;
}