#include "lbkes_k_eigenvalue_solver.h"

#include "B_DO_SteadyState/IterativeOperations/sweep_wgs_context.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_linear_solver.h"

/**Initializes Within-GroupSet solvers.*/
void lbs::DiscOrdKEigenvalueSolver::InitializeWGSSolvers()
{
  wgs_solvers_.clear(); //this is required
  for (auto& groupset : groupsets_)
  {
    std::shared_ptr<SweepChunk> sweep_chunk = SetSweepChunk(groupset);

    auto sweep_wgs_context_ptr =
    std::make_shared<SweepWGSContext<Mat, Vec, KSP>>(
      *this, groupset,
        active_set_source_function_,
        APPLY_WGS_SCATTER_SOURCES,  //lhs_scope
        APPLY_AGS_SCATTER_SOURCES,  //rhs_scope
        true/*with_delayed_psi*/,
        options_.verbose_inner_iterations,
        sweep_chunk);

    auto wgs_solver =
      std::make_shared<WGSLinearSolver<Mat,Vec,KSP>>(sweep_wgs_context_ptr);

    wgs_solvers_.push_back(wgs_solver);
  }//for groupset

}