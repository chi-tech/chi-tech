#include "diffusion_mip.h"

#include "ChiMath/PETScUtils/petsc_utils.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"

#include "ChiPhysics/chi_physics_namespace.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Solves the system and stores the local solution in the vector provide.
 *
 * \param solution Vector in to which the solution will be parsed.*/
void lbs::acceleration::DiffusionMIPSolver::Solve(std::vector<double>& solution)
{
  const std::string fname = "lbs::acceleration::DiffusionMIPSolver::Solve";
  Vec x;
  VecDuplicate(rhs_, &x);
  VecSet(x,0.0);
  KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE);

  KSPSetTolerances(ksp_, 1.e-50,
                   options.residual_tolerance, 1.0e50,
                   options.max_iters);

  if (options.perform_symmetry_check)
  {
    PetscBool symmetry = PETSC_FALSE;
    MatIsSymmetric(A_, 1.0e-6, &symmetry);
    if (symmetry == PETSC_FALSE)
      throw std::logic_error(fname + ":Symmetry check failed");
  }

  if (options.verbose)
  {
    using namespace chi_math::PETScUtils;
    KSPSetConvergenceTest(ksp_, &RelativeResidualConvergenceTest,
                          nullptr, nullptr);

    KSPMonitorSet(ksp_, &KSPMonitorRelativeToRHS, nullptr, nullptr);

    double rhs_norm;
    VecNorm(rhs_, NORM_2, &rhs_norm);
    chi::log.Log() << "RHS-norm " << rhs_norm;
  }


  //============================================= Solve
  KSPSolve(ksp_, rhs_, x);


  //============================================= Print convergence info
  if (options.verbose)
  {
    double sol_norm;
    VecNorm(x, NORM_2, &sol_norm);
    chi::log.Log() << "Solution-norm " << sol_norm;

    using namespace chi_physics;
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp_, &reason);

    chi::log.Log() << "Convergence Reason: "
                   << GetPETScConvergedReasonstring(reason);
  }

  //============================================= Transfer petsc solution to
  //                                              vector
  sdm_.LocalizePETScVector(x, solution, uk_man_);

  //============================================= Cleanup x
  VecDestroy(&x);
}

//###################################################################
/**Solves the system and stores the local solution in the vector provide.
 *
 * \param petsc_solution Vector in to which the solution will be parsed.*/
void lbs::acceleration::DiffusionMIPSolver::Solve(Vec petsc_solution)
{
  const std::string fname = "lbs::acceleration::DiffusionMIPSolver::Solve";
  Vec x;
  VecDuplicate(rhs_, &x);
  VecSet(x,0.0);
  KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE);

  KSPSetTolerances(ksp_, 1.e-50,
                   options.residual_tolerance, 1.0e50,
                   options.max_iters);

  if (options.perform_symmetry_check)
  {
    PetscBool symmetry = PETSC_FALSE;
    MatIsSymmetric(A_, 1.0e-6, &symmetry);
    if (symmetry == PETSC_FALSE)
      throw std::logic_error(fname + ":Symmetry check failed");
  }

  if (options.verbose)
  {
    using namespace chi_math::PETScUtils;
    KSPSetConvergenceTest(ksp_, &RelativeResidualConvergenceTest,
                          nullptr, nullptr);

    KSPMonitorSet(ksp_, &KSPMonitorRelativeToRHS, nullptr, nullptr);

    double rhs_norm;
    VecNorm(rhs_, NORM_2, &rhs_norm);
    chi::log.Log() << "RHS-norm " << rhs_norm;
  }


  //============================================= Solve
  KSPSolve(ksp_, rhs_, x);


  //============================================= Print convergence info
  if (options.verbose)
  {
    double sol_norm;
    VecNorm(x, NORM_2, &sol_norm);
    chi::log.Log() << "Solution-norm " << sol_norm;

    using namespace chi_physics;
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp_, &reason);

    chi::log.Log() << "Convergence Reason: "
                   << GetPETScConvergedReasonstring(reason);
  }

  //============================================= Transfer petsc solution to
  //                                              vector
  VecCopy(x, petsc_solution);

  //============================================= Cleanup x
  VecDestroy(&x);
}