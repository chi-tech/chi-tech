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
  VecDuplicate(m_rhs, &x);
  VecSet(x,0.0);
  KSPSetInitialGuessNonzero(m_ksp,PETSC_FALSE);

  KSPSetTolerances(m_ksp,1.e-50,
                   options.residual_tolerance,1.0e50,
                   options.max_iters);

  if (options.perform_symmetry_check)
  {
    PetscBool symmetry = PETSC_FALSE;
    MatIsSymmetric(m_A, 1.0e-6, &symmetry);
    if (symmetry == PETSC_FALSE)
      throw std::logic_error(fname + ":Symmetry check failed");
  }

  if (options.verbose)
  {
    using namespace chi_math::PETScUtils;
    KSPSetConvergenceTest(m_ksp, &RelativeResidualConvergenceTest,
                          nullptr, nullptr);

    KSPMonitorSet(m_ksp, &GeneralKSPMonitor, nullptr, nullptr);

    double rhs_norm;
    VecNorm(m_rhs, NORM_2, &rhs_norm);
    chi::log.Log() << "RHS-norm " << rhs_norm;
  }


  //============================================= Solve
  KSPSolve(m_ksp,m_rhs,x);


  //============================================= Print convergence info
  if (options.verbose)
  {
    double sol_norm;
    VecNorm(x, NORM_2, &sol_norm);
    chi::log.Log() << "Solution-norm " << sol_norm;

    using namespace chi_physics;
    KSPConvergedReason reason;
    KSPGetConvergedReason(m_ksp, &reason);

    chi::log.Log() << "Convergence Reason: "
                   << GetPETScConvergedReasonstring(reason);
  }

  //============================================= Transfer petsc solution to
  //                                              vector
  m_sdm.LocalizePETScVector(x, solution, m_uk_man);

  //============================================= Cleanup x
  VecDestroy(&x);
}