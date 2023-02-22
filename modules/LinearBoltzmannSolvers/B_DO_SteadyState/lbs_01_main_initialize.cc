#include "lbs_DO_steady_state.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"
#include "IterativeOperations/sweep_wgs_context.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_linear_solver.h"

#include "ChiMPI/chi_mpi.h"
#include "chi_log.h"

//###################################################################
/** Initialize the solver.*/
void lbs::DiscOrdSteadyStateSolver::Initialize()
{
  LBSSolver::Initialize();

  // Initialize source func
  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&LBSSolver::SetSource, this, _1, _2, _3, _4);

  //================================================== Initialize groupsets for
  //                                                   sweeping
  InitializeSweepDataStructures();
  for (auto& groupset : groupsets_)
  {
    InitFluxDataStructures(groupset);

    InitWGDSA(groupset);
    InitTGDSA(groupset);
  }

  InitializeSolverSchemes();           //j
  source_event_tag_ = chi::log.GetRepeatingEventTag("Set Source");
}

/**Initializes Within-GroupSet solvers.*/
void lbs::DiscOrdSteadyStateSolver::InitializeWGSSolvers()
{
  wgs_solvers_.clear(); //this is required
  for (auto& groupset : groupsets_)
  {
    std::shared_ptr<SweepChunk> sweep_chunk = SetSweepChunk(groupset);

    auto sweep_wgs_context_ptr =
    std::make_shared<SweepWGSContext<Mat, Vec, KSP>>(
      *this, groupset,
        active_set_source_function_,
        APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES,  //lhs_scope
        APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES |
        APPLY_AGS_FISSION_SOURCES,                              //rhs_scope
        true/*with_delayed_psi*/,
        options_.verbose_inner_iterations,
        sweep_chunk);

    auto wgs_solver =
      std::make_shared<WGSLinearSolver<Mat,Vec,KSP>>(sweep_wgs_context_ptr);

    wgs_solvers_.push_back(wgs_solver);
  }//for groupset

}