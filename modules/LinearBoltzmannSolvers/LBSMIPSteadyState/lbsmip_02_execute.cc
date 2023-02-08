#include "lbsmip_steady_solver.h"

#include "IterativeOperations/mip_wgs_context.h"
#include "LBSSteadyState/IterativeMethods/wgs_linear_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "ChiConsole/chi_console.h"

#include <iomanip>


/**Execute function.*/
void lbs::MIPSteadyStateSolver::Execute()
{
  MPI_Barrier(MPI_COMM_WORLD);
  for (auto& groupset : groupsets_)
  {
    chi::log.Log()
      << "\n********* Initializing Groupset " << groupset.id
      << "\n" << std::endl;

    InitTGDSA(groupset);

    SolveGroupset(groupset);

    CleanUpTGDSA(groupset);

    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (options_.use_precursors)
    ComputePrecursors();

  UpdateFieldFunctions();

  chi::log.Log() << "LB solver " << TextName() << " execution completed\n";
}

//###################################################################
/**Solves a single groupset.*/
void lbs::MIPSteadyStateSolver::SolveGroupset(LBSGroupset& groupset)
{
  source_event_tag_ = chi::log.GetRepeatingEventTag("Set Source");

  q_moments_local_.assign(q_moments_local_.size(), 0.0);

  auto mip_wgs_context_ptr =
    std::make_shared<MIPWGSContext<Mat, Vec, KSP>>(
      *this, groupset,
      active_set_source_function_,
      APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES |
      SUPPRESS_WG_SCATTER,                                    //lhs_scope
      APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES |
      APPLY_AGS_FISSION_SOURCES,                             //rhs_scope
      options_.verbose_inner_iterations);

  WGSLinearSolver<Mat,Vec,KSP> solver(mip_wgs_context_ptr);
  solver.Setup();
  solver.Solve();

  if (options_.write_restart_data)
    WriteRestartData(options_.write_restart_folder_name,
                     options_.write_restart_file_base);

  chi::log.Log()
    << "Groupset solve complete.                  Process memory = "
    << std::setprecision(3)
    << chi_objects::ChiConsole::GetMemoryUsageInMB() << " MB";
}