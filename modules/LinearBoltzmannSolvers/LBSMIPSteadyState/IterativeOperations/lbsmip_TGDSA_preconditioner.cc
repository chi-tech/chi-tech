#include "lbsmip_shell_operations.h"
#include "LBSMIPSteadyState/Tools/lbsmip_ksp_data_context.h"

#include "LBSSteadyState/Acceleration/diffusion_mip.h"

#include "LBSSteadyState/Groupset/lbs_groupset.h"

//###################################################################
/**Applies TGDSA to the given input vector.*/
int lbs::MIP_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output)
{
  constexpr bool NO_DELAYED_PSI = false;
  MIPKSPDataContext* context;
  PCShellGetContext(pc,&context);

  //Shorten some names
  lbs::MIPSteadyStateSolver& solver = context->solver;
  LBSGroupset& groupset  = context->groupset;
  auto& petsc_phi_delta = context->operating_vector;

  //============================================= Copy PETSc vector to STL
  auto& phi_delta = context->phi_new_local;
  solver.SetPrimarySTLvectorFromGSPETScVec(groupset, phi_input, phi_delta,
                                           NO_DELAYED_PSI);

  //============================================= Apply TGDSA
  if (groupset.apply_tgdsa)
  {
    std::vector<double> delta_phi_local;
    solver.AssembleTGDSADeltaPhiVector(groupset, phi_delta, delta_phi_local);

    groupset.tgdsa_solver->Assemble_b(delta_phi_local);
    groupset.tgdsa_solver->Solve(delta_phi_local);

    solver.DisAssembleTGDSADeltaPhiVector(groupset, delta_phi_local, phi_delta);
  }

  //============================================= Copy STL vector to PETSc Vec
  // We copy the STL data to the operating vector
  // petsc_phi_delta first because its already sized.
  // pc_output is not necessarily initialized yet.
  solver.SetGSPETScVecFromPrimarySTLvector(groupset, petsc_phi_delta, phi_delta,
                                           NO_DELAYED_PSI);

  //============================================= Return result
  VecCopy(petsc_phi_delta, pc_output);

  return 0;
}