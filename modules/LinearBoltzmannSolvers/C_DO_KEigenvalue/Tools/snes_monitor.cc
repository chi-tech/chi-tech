#include "../lbkes_k_eigenvalue_solver.h"
#include "A_LBSSolver/IterativeMethods/ags_context.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include <petscsnes.h>
#include <iomanip>

namespace lbs
{

PetscErrorCode KEigenSNESMonitor(SNES snes, PetscInt iter, PetscReal rnorm, void*)
{
  lbs::AGSContext<Mat,Vec,KSP>* context;
  SNESGetApplicationContext(snes, &context);

  auto& lbs_solver = context->lbs_solver_;

  auto keigen_solver = dynamic_cast<lbs::DiscOrdKEigenvalueSolver*>(&lbs_solver);
  if (not keigen_solver)
    throw std::logic_error("lbs::KEigenSNESMonitor: Solver cast failure");

  std::stringstream iter_info;
  iter_info
    << chi::program_timer.GetTimeString() << " "
    << lbs_solver.TextName()
    << "_NonLinearK_Outer"
    << " Iteration " << std::setw(5) << iter
    << " Residual " << std::setw(11) << rnorm
    << " k_eff " << std::fixed << std::setw(9) << std::setprecision(7)
    << keigen_solver->GetKeff();

  chi::log.Log() << iter_info.str();

  return 0;
}

}//namespace lbs