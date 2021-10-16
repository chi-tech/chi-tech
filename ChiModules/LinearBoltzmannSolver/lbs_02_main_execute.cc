#include "lbs_linear_boltzmann_solver.h"
#include "IterativeMethods/lbs_iterativemethods.h"

#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include "chi_log.h"
extern ChiLog&     chi_log;

#include "chi_mpi.h"
extern ChiMPI&      chi_mpi;

#include "ChiConsole/chi_console.h"
extern ChiConsole&  chi_console;

#include <iomanip>

//###################################################################
/**Execute the solver.*/
void LinearBoltzmann::Solver::Execute()
{
  MPI_Barrier(MPI_COMM_WORLD);
  for (auto& groupset : groupsets)
  {
    chi_log.Log(LOG_0)
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

  if (options.use_precursors)
    ComputePrecursors();

  chi_log.Log(LOG_0) << "NPTransport solver execution completed\n";
}


//###################################################################
/**Solves a single groupset.*/
void LinearBoltzmann::Solver::SolveGroupset(LBSGroupset& groupset)
{
  source_event_tag = chi_log.GetRepeatingEventTag("Set Source");

  //================================================== Setting up required
  //                                                   sweep chunks
  auto sweep_chunk = SetSweepChunk(groupset);
  MainSweepScheduler sweep_scheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
                                     groupset.angle_agg,
                                     *sweep_chunk);

  q_moments_local.assign(q_moments_local.size(), 0.0);

  if (groupset.iterative_method == IterativeMethod::CLASSICRICHARDSON)
  {
    ClassicRichardson(groupset, sweep_scheduler,
                      APPLY_MATERIAL_SOURCE |
                      APPLY_AGS_SCATTER_SOURCE | APPLY_WGS_SCATTER_SOURCE |
                      APPLY_AGS_FISSION_SOURCE | APPLY_WGS_FISSION_SOURCE,
                      options.verbose_inner_iterations);
  }
  else if (groupset.iterative_method == IterativeMethod::GMRES)
  {
    GMRES(groupset, sweep_scheduler,
          APPLY_WGS_SCATTER_SOURCE | APPLY_WGS_FISSION_SOURCE,  //lhs_scope
          APPLY_MATERIAL_SOURCE | APPLY_AGS_SCATTER_SOURCE |
          APPLY_AGS_FISSION_SOURCE,                             //rhs_scope
          options.verbose_inner_iterations);
  }

  if (options.write_restart_data)
    WriteRestartData(options.write_restart_folder_name,
                     options.write_restart_file_base);

  chi_log.Log(LOG_0)
    << "Groupset solve complete.                  Process memory = "
    << std::setprecision(3)
    << chi_console.GetMemoryUsageInMB() << " MB";
}

