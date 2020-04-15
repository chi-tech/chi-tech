#include <ChiLog/chi_log.h>
#include "petsc_utils.h"

#include "chi_log.h"
extern ChiLog chi_log;

//###################################################################
/**Creates a general square matrix.*/
Mat chi_math::PETScUtils::CreateSquareMatrix(int local_size, int global_size)
{
  Mat A;
  MatCreate(PETSC_COMM_WORLD,&A);
  MatSetSizes(A,local_size, local_size,
                global_size, global_size);
  MatSetType(A,MATMPIAIJ);

  return A;
}

//###################################################################
/**Initializes the sparsity pattern of a matrix.*/
void chi_math::PETScUtils::InitMatrixSparsity(
  Mat A,
  std::vector<int>& nodal_nnz_in_diag,
  std::vector<int>& nodal_nnz_off_diag)
{
  MatMPIAIJSetPreallocation(A,0,nodal_nnz_in_diag.data(),
                            0,nodal_nnz_off_diag.data());
  MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
  MatSetUp(A);
}

//###################################################################
/**Creates a general vector.*/
Vec chi_math::PETScUtils::CreateVector(int local_size, int global_size)
{
  Vec x;
  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetSizes(x, local_size, global_size);
  VecSetType(x,VECMPI);

  return x;
}

//###################################################################
/**Creates a common Krylov-solver setup.*/
chi_math::PETScUtils::PETScSolverSetup
  chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
    Mat ref_matrix,
    std::string in_solver_name,
    std::string in_solver_type,
    std::string in_preconditioner_type,
    double in_relative_residual_tolerance,
    int in_maximum_iterations)
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

  KSPMonitorSet(setup.ksp,&GeneralKSPMonitor,NULL,NULL);
  KSPSetConvergenceTest(setup.ksp,&RelativeResidualConvergenceTest,NULL,NULL);

  return setup;
}

//###################################################################
/**Relative Residual Convergence test.*/
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
  if (rhs_norm < 1.0e-25)
    rhs_norm = 1.0;

  //============================== Compute test criterion
  double tol;
  int    maxIts;
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
  if (rhs_norm < 1.0e-25)
    rhs_norm = 1.0;

  char buff[100];
//  char ksp_name[50] = "Nothing";
  const char* ksp_name;
  KSPGetOptionsPrefix(ksp,&ksp_name);
  if (rnorm/rhs_norm < 1.0e-2)
  {
    snprintf(buff,100,"%s iteration %4d - Residual %.3e\n",
      ksp_name,n,rnorm/rhs_norm);
  }
  else
  {
    snprintf(buff,100,"%s iteration %4d - Residual %.7f\n",
      ksp_name,n,rnorm/rhs_norm);
  }

  chi_log.Log(LOG_0) << buff;



  return 0;
}

