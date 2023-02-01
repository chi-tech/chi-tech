#include "lbs_shell_operations.h"
#include "LBSSteadyState/Tools/ksp_data_context.h"

#include "LBSSteadyState/Acceleration/diffusion_mip.h"

#include "LBSSteadyState/Groupset/lbs_groupset.h"

//###################################################################
/**Applies WGDSA or TGDSA to the given input vector.*/
int lbs::WGDSA_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output)
{
  constexpr bool WITH_DELAYED_PSI = true;
  KSPDataContext* context;
  PCShellGetContext(pc,&context);

  //Shorten some names
  lbs::SteadyStateSolver& solver = context->solver;
  LBSGroupset& groupset  = context->groupset;
  auto& petsc_phi_delta = context->operating_vector;

  //============================================= Copy PETSc vector to STL
  auto& phi_delta = solver.phi_new_local;
  solver.SetSTLvectorFromPETScVec(groupset, phi_input, phi_delta,
                                  WITH_DELAYED_PSI);
  groupset.angle_agg.SetDelayedPsiOld2New();

  //============================================= Apply WGDSA
  if (groupset.apply_wgdsa)
  {
    solver.AssembleWGDSADeltaPhiVector(groupset, phi_delta);

    groupset.wgdsa_solver->Assemble_b(solver.delta_phi_local);
    groupset.wgdsa_solver->Solve(solver.delta_phi_local);

    solver.DisAssembleWGDSADeltaPhiVector(groupset, phi_delta);
  }
  //============================================= Apply TGDSA
  if (groupset.apply_tgdsa)
  {
    solver.AssembleTGDSADeltaPhiVector(groupset, phi_delta);

    groupset.tgdsa_solver->Assemble_b(solver.delta_phi_local);
    groupset.tgdsa_solver->Solve(solver.delta_phi_local);

    solver.DisAssembleTGDSADeltaPhiVector(groupset, phi_delta);
  }

  //============================================= Copy STL vector to PETSc Vec
  // We copy the STL data to the operating vector
  // petsc_phi_delta first because its already sized.
  // pc_output is not necessarily initialized yet.
  solver.SetPETScVecFromSTLvector(groupset, petsc_phi_delta, phi_delta,
                                  WITH_DELAYED_PSI);

  //============================================= Return result
  VecCopy(petsc_phi_delta, pc_output);

  return 0;
}