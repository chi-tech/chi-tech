#include "lbs_linear_boltzmann_solver.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"
#include "IterativeOperations/sweep_wgs_context.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_linear_solver.h"

#include "ChiMPI/chi_mpi.h"
#include "chi_log.h"

//###################################################################
/** Initialize the solver.*/
void lbs::SteadyStateSolver::Initialize()
{
  LBSSolver::Initialize();

  // Initialize source func
  using namespace std::placeholders;
  active_set_source_function_ =
    std::bind(&SteadyStateSolver::SetSource, this, _1, _2, _3, _4);

  //================================================== Initialize groupsets_ for
  //                                                   sweeping
  chi::log.Log() << "Initializing groupset sweeping data" << TextName() << ".";
  for (auto& groupset : groupsets_)
  {
    ComputeSweepOrderings(groupset);
    InitFluxDataStructures(groupset);

    InitWGDSA(groupset);
    InitTGDSA(groupset);
  }

  InitializeSolverSchemes();           //j
  source_event_tag_ = chi::log.GetRepeatingEventTag("Set Source");
}

void lbs::SteadyStateSolver::InitializeSolverSchemes()
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

/**Initializes Within-GroupSet solvers.*/
void lbs::SteadyStateSolver::InitializeWGSSolvers()
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