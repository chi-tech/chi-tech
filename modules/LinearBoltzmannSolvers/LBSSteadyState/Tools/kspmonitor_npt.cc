#include "kspmonitor_npt.h"
#include "LBSSteadyState/lbs_linear_boltzmann_solver.h"

#include "ksp_data_context.h"

#include <iomanip>
#include "ChiLog/chi_log.h"
#include "ChiTimer/chi_timer.h"

//###################################################################
/**Customized convergence test.*/
PetscErrorCode lbs::
  KSPConvergenceTest(KSP ksp, PetscInt n, PetscReal rnorm,
                        KSPConvergedReason* convergedReason, void*)
{
  constexpr bool WITH_DELAYED_PSI = true;
  //======================================== Get data context
  KSPDataContext* context;
  KSPGetApplicationContext(ksp,&context);

  //======================================== Set rhs norm
  auto rhs_norm = context->rhs_preconditioned_norm;
  if (rhs_norm < 1.0e-25)
    rhs_norm = 1.0;

  //======================================== Compute test criterion
  double tol;
  int64_t    maxIts;
  KSPGetTolerances(ksp, nullptr,&tol, nullptr,&maxIts);

  double relative_residual = rnorm/rhs_norm;

  //======================================== Print iteration information
  std::string offset;
  if (context->groupset.apply_wgdsa || context->groupset.apply_tgdsa)
    offset = std::string("    ");

  std::stringstream iter_info;
  iter_info
    << chi::program_timer.GetTimeString() << " "
    << offset
    << "WGS groups ["
    << context->groupset.groups.front().id
    << "-"
    << context->groupset.groups.back().id
    << "]"
    << " Iteration " << std::setw(5) << n
    << " Residual " << std::setw(9) << relative_residual;

  if (relative_residual < tol)
  {
    *convergedReason = KSP_CONVERGED_RTOL;
    iter_info << " CONVERGED\n";
  }

  if (context->solver.Options().verbose_inner_iterations)
    chi::log.Log() << iter_info.str() << std::endl;

  const double SIXTY_SECOND_INTERVAL = 60000.0; //time in milliseconds
  if (context->groupset.iterative_method == IterativeMethod::KRYLOV_GMRES or
      context->groupset.iterative_method == IterativeMethod::KRYLOV_GMRES_CYCLES)
  {
    if (context->last_iteration == n)
    {
      if (context->solver.Options().write_restart_data)
      {
        if ((chi::program_timer.GetTime()/SIXTY_SECOND_INTERVAL) >
          context->solver.LastRestartWrite() +
          context->solver.Options().write_restart_interval)
        {
          Vec phi_new;
          KSPBuildSolution(ksp, nullptr,&phi_new);

          context->solver.
            SetPrimarySTLvectorFromGSPETScVec(context->groupset, phi_new,
                                              context->phi_old_local,
                                              WITH_DELAYED_PSI);

          context->solver.LastRestartWrite() =
            chi::program_timer.GetTime()/SIXTY_SECOND_INTERVAL;
          context->solver.WriteRestartData(
            context->solver.Options().write_restart_folder_name,
            context->solver.Options().write_restart_file_base);
        }
      }
    }
  }
  context->last_iteration = n;

  return KSP_CONVERGED_ITERATING;
}
