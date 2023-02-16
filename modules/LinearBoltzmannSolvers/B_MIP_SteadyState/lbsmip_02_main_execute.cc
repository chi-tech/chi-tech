#include "lbsmip_steady_solver.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"

#include "chi_log.h"

//###################################################################
/**Execute the solver.*/
void lbs::MIPSteadyStateSolver::Execute()
{
  primary_ags_solver_->Setup();
  primary_ags_solver_->Solve();

  if (options_.use_precursors)
    ComputePrecursors();

  UpdateFieldFunctions();

  chi::log.Log() << "LB solver " << TextName() << " execution completed\n";
}
