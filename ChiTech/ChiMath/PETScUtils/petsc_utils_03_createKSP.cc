#include "petsc_utils.h"

#include <iomanip>

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Creates a common Krylov-solver setup.
 *
This is a macro for:
\code
PETScSolverSetup setup;

KSPCreate(PETSC_COMM_WORLD,&setup.ksp);
KSPSetOperators(setup.ksp,ref_matrix,ref_matrix);
KSPSetType(setup.ksp,in_solver_type.c_str());

KSPSetOptionsPrefix(setup.ksp,in_solver_name.c_str());

KSPGetPC(setup.ksp,&setup.pc);
PCSetType(setup.pc,in_preconditioner_type.c_str());

KSPSetTolerances(setup.ksp,1.e-50,
                 in_relative_residual_tolerance,1.0e50,
                 in_maximum_iterations);
KSPSetInitialGuessNonzero(setup.ksp,PETSC_TRUE);

KSPMonitorSet(setup.ksp,&GeneralKSPMonitor,NULL,NULL);
KSPSetConvergenceTest(setup.ksp,&RelativeResidualConvergenceTest,NULL,NULL);

return setup;
\endcode*/
chi_math::PETScUtils::PETScSolverSetup
chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
  Mat ref_matrix,
  std::string in_solver_name,
  std::string in_solver_type,
  std::string in_preconditioner_type,
  double in_relative_residual_tolerance,
  int64_t in_maximum_iterations)
{
  PETScSolverSetup setup;

  KSPCreate(PETSC_COMM_WORLD,&setup.ksp);
  KSPSetOperators(setup.ksp,ref_matrix,ref_matrix);
  KSPSetType(setup.ksp,in_solver_type.c_str());

  KSPSetOptionsPrefix(setup.ksp,in_solver_name.c_str());

  KSPGetPC(setup.ksp,&setup.pc);
  PCSetType(setup.pc,in_preconditioner_type.c_str());

  KSPSetTolerances(setup.ksp,1.e-50,
                   in_relative_residual_tolerance,1.0e50,
                   in_maximum_iterations);
  KSPSetInitialGuessNonzero(setup.ksp,PETSC_TRUE);

//  KSPMonitorSet(setup.ksp,&GeneralKSPMonitor,NULL,NULL);
  KSPSetConvergenceTest(setup.ksp,&RelativeResidualConvergenceTest,NULL,NULL);
  KSPSetFromOptions(setup.ksp);

  KSPMonitorSet(setup.ksp,&chi_math::PETScUtils::GeneralKSPMonitor,NULL,NULL);

  return setup;
}

//###################################################################
/**Relative Residual Convergence test. The test uses
 * the L2-norm of the residual (\f$ ||b-Ax||_2 \f$), divided by
 * by the L2-norm of the right hand side (\f$ ||b||_2 \f$) compared
 * to a tolerance, \f$ \epsilon \f$.

\f[
 \frac{||b-Ax||_2}{||b||_2} < \epsilon
\f]
 *
 * */
PetscErrorCode
chi_math::PETScUtils::RelativeResidualConvergenceTest(
  KSP ksp, PetscInt n,
  PetscReal rnorm,
  KSPConvergedReason* convergedReason,
  void *monitordestroy)
{
  //============================== Compute rhs norm
  Vec Rhs;
  KSPGetRhs(ksp,&Rhs);
  double rhs_norm;
  VecNorm(Rhs,NORM_2,&rhs_norm);
  if (rhs_norm < 1.0e-12)
    rhs_norm = 1.0;

  //============================== Compute test criterion
  double tol;
  int64_t maxIts;
  KSPGetTolerances(ksp,NULL,&tol,NULL,&maxIts);

  double relative_residual = rnorm/rhs_norm;

  if (relative_residual < tol)
    *convergedReason = KSP_CONVERGED_RTOL;

  return KSP_CONVERGED_ITERATING;
}

//###################################################################
/**General monitor.*/
PetscErrorCode chi_math::PETScUtils::GeneralKSPMonitor(
  KSP ksp, PetscInt n, PetscReal rnorm, void *monitordestroy)
{
  Vec Rhs;
  KSPGetRhs(ksp,&Rhs);
  double rhs_norm;
  VecNorm(Rhs,NORM_2,&rhs_norm);
  if (rhs_norm < 1.0e-12)
    rhs_norm = 1.0;

  const char* ksp_name;
  KSPGetOptionsPrefix(ksp,&ksp_name);

  std::stringstream buff;
  buff
    << ksp_name
    << " iteration "
    << std::setw(4) << n
    << " - Residual "
    << std::scientific << std::setprecision(7) << rnorm / rhs_norm
    << std::endl;

  chi_log.Log(LOG_0) << buff.str();

  return 0;
}