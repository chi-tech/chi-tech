#include "A_LBSSolver/lbs_solver.h"
#include "ags_linear_solver.h"
#include "wgs_context.h"
#include "snes_k_monitor.h"
#include "snes_k_residual_func_context.h"
#include "snes_k_residual_function.h"

#include "ChiMath/PETScUtils/petsc_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <iomanip>

#include <petscsnes.h>

#define sc_int64_t static_cast<int64_t>

namespace lbs
{

/**This routine applies PETSc's Scalabel Non-linear Equation Solver (SNES) to
 * to solve the k-eigenvalue problem using PJFNK.*/
int NonLinearKEigen(LBSSolver& lbs_solver,
                    double nonlinear_abs_tolerance,
                    int    nonlinear_max_iterations,
                    double& k_eff)
{
  chi::log.Log()
    << "\n********** Solving k-eigenvalue problem with "
    << "Non-Linear PJFNK.\n";

  std::vector<int> groupset_ids;
  for (auto& groupset : lbs_solver.Groupsets())
    groupset_ids.push_back(groupset.id_);

  const auto [num_local_dofs, num_globl_dofs] =
    lbs_solver.GetNumPhiIterativeUnknowns();

  //============================================= Create a context to keep
  //                                              track of k_eff
  const std::string solver_name = lbs_solver.TextName();
  KResidualFunctionContext residual_function_context{solver_name, /*k_eff=*/1.0};

  //============================================= Create the vectors
  Vec phi, r; /* solution, residual vectors */
  phi = chi_math::PETScUtils::CreateVector(sc_int64_t(num_local_dofs),
                                           sc_int64_t(num_globl_dofs));
  VecDuplicate(phi, &r);

  auto primary_ags_solver = lbs_solver.GetPrimaryAGSSolver();

  //============================================= Create SNES
  SNES        snes; /* nonlinear solver context */

  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESSetType(snes, SNESNEWTONLS);
  SNESSetApplicationContext(snes, &(*primary_ags_solver->GetContext()));
//  SNESSetFromOptions(snes);

  SNESSetTolerances(snes,
                    nonlinear_abs_tolerance,  //absolute tolerance
                    1.0e-50,                  //relative tolerance
                    nonlinear_abs_tolerance,  //solution tolerance
                    nonlinear_max_iterations, //max iterations
                    -1);                      //max r-evals

  SNESSetMaxLinearSolveFailures(snes, nonlinear_max_iterations);

  SNESMonitorSet(snes, &lbs::KEigenSNESMonitor,
                 &residual_function_context, nullptr);

  //============================================= Set the residual function
  SNESSetFunction(snes, r, SNESKResidualFunction, &residual_function_context);

  //============================================= Setup Matrix-Free Jacobian
  Mat         J;    /* Jacobian matrix */

  MatCreateSNESMF(snes, &J);
  SNESSetJacobian(snes, J, J, MatMFFDComputeJacobian, nullptr);

  //============================================= Set linear solver and
  //                                              preconditioner
  KSP         ksp;  /* linear solver context */
  PC          pc;   /* preconditioner context */

  SNESGetKSP(snes, &ksp);
  KSPSetType(ksp, KSPGMRES);

  auto& first_groupset = lbs_solver.Groupsets().front();
  KSPSetTolerances(ksp,
                   first_groupset.residual_tolerance_, //relative tol
                   first_groupset.residual_tolerance_, //absolute tol
                   1.0e6,                              //divergence tol
                   first_groupset.max_iterations_);    //max iters


  KSPGMRESSetRestart(ksp, first_groupset.gmres_restart_intvl_);

  KSPSetOptionsPrefix(ksp, ("    " + solver_name + "_NonLinearK_Inner").c_str());
  if (lbs_solver.Options().verbose_inner_iterations)
    KSPMonitorSet(ksp, &chi_math::PETScUtils::KSPMonitorStraight, nullptr, nullptr);

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCNONE);

  //============================================= Compute initial guess
  lbs_solver.SetPhiVectorScalarValues(lbs_solver.PhiOldLocal(), 1.0);

  lbs_solver.SetMultiGSPETScVecFromPrimarySTLvector(
    groupset_ids, phi, PhiSTLOption::PHI_OLD);

  //============================================= Solve
  SNESSolve(snes, nullptr, phi);

  //============================================= Unpack solution
  const auto& groups = lbs_solver.Groups();
  lbs_solver.SetPrimarySTLvectorFromGroupScopedPETScVec(
    groups.front().id_, groups.back().id_, phi, lbs_solver.PhiOldLocal());

  //============================================= Compute final k_eff
  k_eff = lbs_solver.ComputeFissionProduction(lbs_solver.PhiOldLocal());

  PetscInt number_of_func_evals;
  SNESGetNumberFunctionEvals(snes, &number_of_func_evals);

  //================================================== Print summary
  chi::log.Log() << "\n"
    << "        Final k-eigenvalue    :        "
    << std::fixed << std::setw(10) << std::setprecision(7)
    << k_eff
    << " (" << number_of_func_evals << ")"
    << "\n";

  return 0;
}

}//namespace lbs