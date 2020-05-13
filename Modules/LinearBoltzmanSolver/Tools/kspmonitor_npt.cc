#include "kspmonitor_npt.h"
#include "../lbs_linear_boltzman_solver.h"

#include "ksp_data_context.h"

#include <iomanip>
#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

extern ChiLog&     chi_log;
extern ChiTimer   chi_program_timer;

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
  KSPDataContext* context;
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

  context->groupset->latest_convergence_metric = std::min(relative_residual, 1.0);

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

  if (context->groupset->iterative_method == NPT_GMRES)
  {
    if (context->last_iteration == n)
    {
      if (context->solver->options.write_restart_data)
      {
        if ((chi_program_timer.GetTime()/60000.0) >
          context->solver->last_restart_write +
          context->solver->options.write_restart_interval)
        {
          Vec phi_new;
          KSPBuildSolution(ksp,NULL,&phi_new);

          context->solver->
          DisAssembleVector(context->groupset, phi_new,
                            context->solver->phi_old_local.data());

          context->solver->last_restart_write = chi_program_timer.GetTime()/60000.0;
          context->solver->WriteRestartData(
            context->solver->options.write_restart_folder_name,
            context->solver->options.write_restart_file_base);
        }
      }
    }
  }
  context->last_iteration = n;

  return KSP_CONVERGED_ITERATING;
}
