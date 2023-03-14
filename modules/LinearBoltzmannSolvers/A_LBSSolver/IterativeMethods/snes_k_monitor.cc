#include "snes_k_residual_func_context.h"
#include "ags_context.h"

#include "A_LBSSolver/lbs_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include <petscsnes.h>
#include <iomanip>

namespace lbs
{

PetscErrorCode KEigenSNESMonitor(SNES snes, PetscInt iter,
                                 PetscReal rnorm, void*ctx)
{
  auto& residual_context = *(KResidualFunctionContext*)ctx;

  double k_eff = residual_context.k_eff;
  double reactivity = (k_eff - 1.0) / k_eff;

  std::stringstream iter_info;
  iter_info
    << chi::program_timer.GetTimeString() << " "
    << residual_context.solver_name
    << "_NonLinearK_Outer"
    << " Iteration " << std::setw(5) << iter
    << " Residual " << std::setw(11) << rnorm
    << " k_eff " << std::fixed << std::setw(10) << std::setprecision(7)
    << k_eff
    << std::setprecision(2)
    << "  reactivity " << std::setw(10) << reactivity * 1e5;

  chi::log.Log() << iter_info.str();

  return 0;
}

}//namespace lbs