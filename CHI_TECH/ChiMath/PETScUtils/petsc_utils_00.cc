#include "petsc_utils.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Creates a general square matrix.
 *
 * This is a macro for:
\code
Mat A;
MatCreate(PETSC_COMM_WORLD,&A);
MatSetType(A,MATMPIAIJ);
MatSetSizes(A,local_size, local_size,
              global_size, global_size);

 return A;
\endcode

*/
Mat chi_math::PETScUtils::CreateSquareMatrix(int local_size, int global_size)
{
  Mat A;
  MatCreate(PETSC_COMM_WORLD,&A);
  MatSetType(A,MATMPIAIJ);
  MatSetSizes(A,local_size, local_size,
                global_size, global_size);

  return A;
}

//###################################################################
/**Creates a general square matrix.
 *
 * This is a macro for:
\code
MatCreate(PETSC_COMM_WORLD,&A);
MatSetType(A,MATMPIAIJ);
MatSetSizes(A,local_size, local_size,
              global_size, global_size);
\endcode

*/
void chi_math::PETScUtils::
  CreateSquareMatrix(Mat& A, int local_size, int global_size)
{
  MatCreate(PETSC_COMM_WORLD,&A);
  MatSetType(A,MATMPIAIJ);
  MatSetSizes(A,local_size, local_size,
              global_size, global_size);
}

//###################################################################
/**Initializes the sparsity pattern of a matrix.

This is a macro for:
\code
MatMPIAIJSetPreallocation(A,0,nodal_nnz_in_diag.data(),
                            0,nodal_nnz_off_diag.data());
MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
MatSetUp(A);
\endcode
*/
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
/**Creates a general vector.
 *
This is a macro for:
\code
Vec x;
VecCreate(PETSC_COMM_WORLD,&x);
VecSetType(x,VECMPI);
VecSetSizes(x, local_size, global_size);
VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

return x;
\endcode*/
Vec chi_math::PETScUtils::
  CreateVector(int local_size, int global_size)
{
  Vec x;
  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetType(x,VECMPI);
  VecSetSizes(x, local_size, global_size);
  VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

  return x;
}

//###################################################################
/**Creates a general vector.
 *
This is a macro for:
\code
VecCreate(PETSC_COMM_WORLD,&x);
VecSetType(x,VECMPI);
VecSetSizes(x, local_size, global_size);
VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
\endcode*/
void chi_math::PETScUtils::
  CreateVector(Vec& x, int local_size, int global_size)
{
  VecCreate(PETSC_COMM_WORLD,&x);
  VecSetType(x,VECMPI);
  VecSetSizes(x, local_size, global_size);
  VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
}

//###################################################################
/**Creates a general vector with ghost value support.
 *
This is a macro for:
\code
Vec x;
VecCreateGhost(PETSC_COMM_WORLD,
               local_size,
               global_size,
               nghosts,
               ghost_indices.data(),
               &x);

VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

return x;
\endcode*/
Vec chi_math::PETScUtils::
  CreateVectorWithGhosts(int local_size, int global_size,
                         int nghosts,
                         std::vector<int>& ghost_indices)
{
  Vec x;
  VecCreateGhost(PETSC_COMM_WORLD,
                 local_size,
                 global_size,
                 nghosts,
                 (ghost_indices.empty())? NULL : ghost_indices.data(),
                 &x);

  VecSetOption(x,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);

  return x;
}

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

//###################################################################
/**Copies a PETSc vector to a STL vector. Only the local portion is
 * copied.*/
void chi_math::PETScUtils::CopyVecToSTLvector(
  Vec x, std::vector<double>& data, size_t N)
{
  data.clear();
  data.resize(N,0.0);
  const double* x_ref;
  VecGetArrayRead(x,&x_ref);

  for (size_t i=0; i<N; ++i)
    data[i] = x_ref[i];

  VecRestoreArrayRead(x,&x_ref);
}

//###################################################################
/**Copies global values from a PETSc vector to a STL vector.*/
void chi_math::PETScUtils::CopyGlobalVecToSTLvector(
  Vec x,
  const std::vector<int>& global_indices,
  std::vector<double> &data)
{
  //=================================== Populating local indices
  size_t N = global_indices.size();
  std::vector<int> local_indices(N,0);
  unsigned int counter=0;
  for (unsigned int val : global_indices)
  {
    local_indices[counter] = counter;
    ++counter;
  }

  //=================================== Creating PETSc vector
  Vec local_vec;
  VecCreateSeq(PETSC_COMM_SELF,global_indices.size()+1,&local_vec);
  VecSet(local_vec,0.0);

  //=================================== Create and transfer index sets
  IS global_set;
  IS local_set;
  ISCreateGeneral(PETSC_COMM_WORLD, N, global_indices.data(),
                  PETSC_COPY_VALUES,&global_set);
  ISCreateGeneral(PETSC_COMM_WORLD, N, local_indices.data(),
                  PETSC_COPY_VALUES,&local_set);
  VecScatter scat;
  VecScatterCreate(x,global_set,local_vec,local_set,&scat);
  VecScatterBegin(scat,x,local_vec,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(scat,x,local_vec,INSERT_VALUES,SCATTER_FORWARD);

  //=================================== Copy to STL
  data.clear();
  data.resize(N,0.0);
  const double* x_ref;
  VecGetArrayRead(x,&x_ref);

  for (size_t i=0; i<N; ++i)
    data[i] = x_ref[i];

  VecRestoreArrayRead(x,&x_ref);

  //=================================== Cleanup
  ISDestroy(&global_set);
  ISDestroy(&local_set);

  VecDestroy(&local_vec);
}

