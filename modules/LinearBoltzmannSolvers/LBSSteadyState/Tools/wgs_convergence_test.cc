#include "wgs_convergence_test.h"

#include "wgs_context.h"
#include "LinearBoltzmannSolvers/LBSSteadyState/Groupset/lbs_groupset.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiTimer/chi_timer.h"

#include <iomanip>

namespace lbs
{

//###################################################################
/**Customized convergence test.*/
PetscErrorCode GSConvergenceTest(KSP ksp, PetscInt n, PetscReal rnorm,
                                 KSPConvergedReason* convergedReason, void*)
{
  //======================================== Get data context
  WGSContext<Mat, Vec, KSP>* context;
  KSPGetApplicationContext(ksp,&context);

  //======================================== Set rhs norm
  double residual_scale = 1.0;
  switch (context->residual_scale_type_)
  {
    case chi_math::ResidualScaleType::NONE:
      residual_scale = 1.0;
      break;
    case chi_math::ResidualScaleType::RHS_NORM:
      if (context->rhs_norm_ > 1.0e-25)
        residual_scale = 1.0/context->rhs_norm_;
      break;
    case chi_math::ResidualScaleType::RHS_PRECONDITIONED_NORM:
      if (context->rhs_preconditioned_norm_ > 1.0e-25)
        residual_scale = 1.0/context->rhs_preconditioned_norm_;
      break;
    case chi_math::ResidualScaleType::CUSTOM_SCALE:
      if (context->custom_residual_scale_ > 1.0e-25)
        residual_scale = 1.0/context->custom_residual_scale_;
      break;
  }

  //======================================== Compute test criterion
  double tol;
  int64_t    maxIts;
  KSPGetTolerances(ksp, nullptr,&tol, nullptr,&maxIts);

  double scaled_residual = rnorm * residual_scale;

  //======================================== Print iteration information
  std::string offset;
  if (context->groupset_.apply_wgdsa || context->groupset_.apply_tgdsa)
    offset = std::string("    ");

  std::stringstream iter_info;
  iter_info
    << chi::program_timer.GetTimeString() << " "
    << offset
    << "WGS groups ["
    << context->groupset_.groups.front().id
    << "-"
    << context->groupset_.groups.back().id
    << "]"
    << " Iteration " << std::setw(5) << n
    << " Residual " << std::setw(9) << scaled_residual;

  if (scaled_residual < tol)
  {
    *convergedReason = KSP_CONVERGED_RTOL;
    iter_info << " CONVERGED\n";
  }

  if (context->log_info_)
    chi::log.Log() << iter_info.str() << std::endl;

  return KSP_CONVERGED_ITERATING;
}

}//namespace lbs