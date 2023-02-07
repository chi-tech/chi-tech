#include "lbsmip_shell_operations.h"
#include "LBSMIPSteadyState/Tools/lbsmip_ksp_data_context.h"
#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include "LBSSteadyState/Acceleration/diffusion_mip.h"

#include "LBSSteadyState/Groupset/lbs_groupset.h"

typedef chi_mesh::sweep_management::SweepScheduler MainSweepScheduler;
//###################################################################
/**Computes the action of the transport matrix on a vector.*/
int lbs::MIPMatrixAction_Ax(Mat matrix, Vec krylov_vector, Vec Av)
{
  constexpr bool NO_DELAYED_PSI = false;
  MIPKSPDataContext* context;
  MatShellGetContext(matrix,&context);

  //Shorten some names
  lbs::MIPSteadyStateSolver& solver = context->solver;
  LBSGroupset& groupset  = context->groupset;

  SourceFlags& lhs_source_scope = context->lhs_scope;
  auto& set_source_function = context->set_source_function;

  //============================================= Copy krylov vector into local
  solver.SetPrimarySTLvectorFromGSPETScVec(groupset,
                                           krylov_vector,
                                           context->phi_old_local, NO_DELAYED_PSI);

  //============================================= Setting the source using
  //                                              updated phi_old
  context->q_moments_local.assign(context->q_moments_local.size(), 0.0);
  set_source_function(groupset, context->q_moments_local,
                      solver.PhiOldLocal(),
                      lhs_source_scope);

  //================================================== Setup GS vectors
  const size_t num_gs_local_dofs = solver.LocalNodeCount()*groupset.groups.size();
  std::vector<double> gs_q_moments_local_(num_gs_local_dofs, 0.0);
  auto gs_phi_new_local_ = gs_q_moments_local_;

  //============================================= Sweeping the new source
  auto& mip_solver = *solver.gs_mip_solvers_[groupset.id];

  solver.SetGSSTLvectorFromPrimarySTLvector(groupset, gs_q_moments_local_,
                                            context->q_moments_local);

  mip_solver.Assemble_b(gs_q_moments_local_);
  mip_solver.Solve(gs_phi_new_local_);

  solver.SetPrimarySTLvectorFromGSSTLvector(groupset, gs_phi_new_local_,
                                            context->phi_new_local);

  //=================================================== Copy local into
  //                                                    operating vector
  // We copy the STL data to the operating vector
  // petsc_phi_delta first because it's already sized.
  // pc_output is not necessarily initialized yet.
  solver.SetGSPETScVecFromPrimarySTLvector(groupset,
                                           context->operating_vector,
                                           context->phi_new_local, NO_DELAYED_PSI);

  //============================================= Computing action
  // A  = [I - DLinvMS]
  // Av = [I - DLinvMS]v
  //    = v - DLinvMSv
  VecWAXPY(Av, -1.0, context->operating_vector, krylov_vector);

  return 0;
}