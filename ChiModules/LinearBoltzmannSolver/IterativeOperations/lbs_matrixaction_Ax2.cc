#include "lbs_shell_operations.h"
#include "LinearBoltzmannSolver/Tools/ksp_data_context.h"
#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include "LinearBoltzmannSolver/Groupset/lbs_groupset.h"

typedef chi_mesh::sweep_management::SweepScheduler MainSweepScheduler;
//###################################################################
/**Computes the action of the transport matrix on a vector.*/
int lbs::MatrixAction_Ax(Mat matrix, Vec krylov_vector, Vec Av)
{
  constexpr bool WITH_DELAYED_PSI = true;
  KSPDataContext* context;
  MatShellGetContext(matrix,&context);

  //Shorten some names
  lbs::SteadySolver& solver = context->solver;
  LBSGroupset& groupset  = context->groupset;
  MainSweepScheduler& sweep_scheduler = context->sweep_scheduler;
  auto& sweep_chunk = context->sweep_scheduler.GetSweepChunk();
  SourceFlags& lhs_source_scope = context->lhs_scope;

  //============================================= Copy krylov vector into local
  solver.SetSTLvectorFromPETScVec(groupset,
                                  krylov_vector,
                                  solver.phi_old_local, WITH_DELAYED_PSI);

  //============================================= Setting the source using
  //                                              updated phi_old
  solver.q_moments_local.assign(solver.q_moments_local.size(), 0.0);
  solver.SetSource(groupset, solver.q_moments_local, lhs_source_scope);

  //============================================= Sweeping the new source
  sweep_chunk.ZeroFluxDataStructures();
  sweep_scheduler.Sweep();

  //=================================================== Copy local into
  //                                                    operating vector
  // We copy the STL data to the operating vector
  // petsc_phi_delta first because its already sized.
  // pc_output is not necessarily initialized yet.
  solver.SetPETScVecFromSTLvector(groupset,
                                  context->operating_vector,
                                  solver.phi_new_local, WITH_DELAYED_PSI);

  //============================================= Computing action
  // A  = [I - DLinvMS]
  // Av = [I - DLinvMS]v
  //    = v - DLinvMSv
  VecWAXPY(Av, -1.0, context->operating_vector, krylov_vector);

  return 0;
}