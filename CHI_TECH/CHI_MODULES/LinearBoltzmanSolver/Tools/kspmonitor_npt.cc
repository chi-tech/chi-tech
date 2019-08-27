#include "kspmonitor_npt.h"
#include "CHI_MODULES/LinearBoltzmanSolver/lbs_linear_boltzman_solver.h"

#include "ksp_data_context.h"

#include <iomanip>
#include <chi_log.h>
#include <CHI_TIMER/chi_timer.h>

extern CHI_LOG     chi_log;
extern CHI_TIMER   chi_program_timer;

//###################################################################
/**Customized monitor for PETSc Krylov sub-space solvers.*/
PetscErrorCode
KSPMonitorNPT(KSP ksp, PetscInt n, PetscReal rnorm, void *monitordestroy)
{
  Vec Rhs;

  KSPGetRhs(ksp,&Rhs);

  double rhs_norm;
  VecNorm(Rhs,NORM_2,&rhs_norm);

  chi_log.Log(LOG_0) << "Iteration " << n << " Residual " << rnorm/rhs_norm;

  return 0;
}



//###################################################################
/**Customized convergence test.*/
PetscErrorCode
KSPConvergenceTestNPT(KSP ksp, PetscInt n, PetscReal rnorm,
                      KSPConvergedReason* convergedReason, void *monitordestroy)
{
  //======================================== Get data context
  KSP_DATA_CONTEXT* context;
  KSPGetApplicationContext(ksp,&context);

  //======================================== Compute rhs norm
  Vec Rhs;
  KSPGetRhs(ksp,&Rhs);
  double rhs_norm;
  VecNorm(Rhs,NORM_2,&rhs_norm);
  if (rhs_norm < 1.0e-25)
    rhs_norm = 1.0;

  //======================================== Compute test criterion
  double tol;
  int    maxIts;
  KSPGetTolerances(ksp,NULL,&tol,NULL,&maxIts);

  double relative_residual = rnorm/rhs_norm;


  //======================================== Print iteration information
  std::string offset;
  if (context->groupset->apply_wgdsa || context->groupset->apply_tgdsa)
    offset = std::string("    ");

  std::stringstream iter_info;
  iter_info
    << chi_program_timer.GetTimeString() << " "
    << offset
    << "WGS groups ["
    << context->groupset->groups.front()->id
    << "-"
    << context->groupset->groups.back()->id
    << "]"
    << " Iteration " << std::setw(5) << n
    << " Residual " << std::setw(9) << relative_residual;

  if (relative_residual < tol)
  {
    *convergedReason = KSP_CONVERGED_RTOL;
    iter_info << " CONVERGED\n";
  }


  chi_log.Log(LOG_0) << iter_info.str() << std::endl;

  return KSP_CONVERGED_ITERATING;
}
