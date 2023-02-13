#include "lbs_solver.h"

#include "A_LBSSolver/Tools/ags_context.h"
#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

void lbs::LBSSolver::InitializeSolverSchemes()
{
  chi::log.Log() << "Initializing Solver schemes";

  InitializeWGSSolvers();

  /*This default behavior covers the situation when no Across-GroupSet (AGS)
   * solvers have been created for this solver.*/
  if (ags_solvers_.empty())
  {
    auto ags_context = std::make_shared<AGSContext<Mat,Vec,KSP>>(
      *this, wgs_solvers_);

    auto ags_solver = std::make_shared<AGSLinearSolver<Mat,Vec,KSP>>(
      "richardson", ags_context,groupsets_.front().id, groupsets_.back().id);
    ags_solver->ToleranceOptions().maximum_iterations_ = 1;

    ags_solvers_.push_back(ags_solver);

    primary_ags_solver_ = ags_solvers_.front();
  }//if ags_solvers.empty()
}