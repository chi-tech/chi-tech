#include "lbs_linear_boltzmann_solver.h"

#include "IterativeOperations/sweep_wgs_context.h"
#include "IterativeMethods/wgs_linear_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_mpi.h"

#include "ChiConsole/chi_console.h"

#include <iomanip>


//###################################################################
/**Execute the solver.*/
void lbs::SteadyStateSolver::Execute()
{
  MPI_Barrier(MPI_COMM_WORLD);
  for (auto& groupset : groupsets_)
  {
    chi::log.Log()
      << "\n********* Initializing Groupset " << groupset.id
      << "\n" << std::endl;

    ComputeSweepOrderings(groupset);
    InitFluxDataStructures(groupset);

    InitWGDSA(groupset);
    InitTGDSA(groupset);

    SolveGroupset(groupset);

    CleanUpWGDSA(groupset);
    CleanUpTGDSA(groupset);

    ResetSweepOrderings(groupset);

    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (options_.use_precursors)
    ComputePrecursors();

  UpdateFieldFunctions();

  chi::log.Log() << "LB solver " << TextName() << " execution completed\n";
}


void lbs::SteadyStateSolver::SolveGroupset(LBSGroupset& groupset)
{
  source_event_tag_ = chi::log.GetRepeatingEventTag("Set Source");

  //================================================== Setting up required
  //                                                   sweep chunks
  q_moments_local_.assign(q_moments_local_.size(), 0.0);
  auto sweep_chunk = SetSweepChunk(groupset);

  if (groupset.iterative_method == IterativeMethod::CLASSICRICHARDSON or
      groupset.iterative_method == IterativeMethod::KRYLOV_RICHARDSON or
      groupset.iterative_method == IterativeMethod::KRYLOV_GMRES or
      groupset.iterative_method == IterativeMethod::KRYLOV_BICGSTAB)
  {
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

    WGSLinearSolver<Mat,Vec,KSP> solver(sweep_wgs_context_ptr);
    solver.Setup();
    solver.Solve();
  }

  if (options_.write_restart_data)
    WriteRestartData(options_.write_restart_folder_name,
                     options_.write_restart_file_base);

  chi::log.Log()
    << "Groupset solve complete.                  Process memory = "
    << std::setprecision(3)
    << chi_objects::ChiConsole::GetMemoryUsageInMB() << " MB";
}
