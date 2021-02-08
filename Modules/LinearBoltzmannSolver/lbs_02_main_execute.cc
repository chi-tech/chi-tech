#include "lbs_linear_boltzmann_solver.h"
#include "IterativeMethods/lbs_iterativemethods.h"

#include "ChiMesh/SweepUtilities/SweepScheduler/sweepscheduler.h"

#include <chi_mpi.h>
#include <chi_log.h>
#include <ChiConsole/chi_console.h>

extern ChiMPI&      chi_mpi;
extern ChiLog&     chi_log;
extern ChiConsole&  chi_console;

#include <iomanip>


//###################################################################
/**Execute the solver.*/
void LinearBoltzmann::Solver::Execute()
{
  MPI_Barrier(MPI_COMM_WORLD);
  int gs=-1;
  for (auto& groupset : group_sets)
  {
    ++gs;
    chi_log.Log(LOG_0)
      << "\n********* Initializing Groupset " << gs << "\n" << std::endl;

    groupset.BuildDiscMomOperator(options.scattering_order,
                                  options.geometry_type);
    groupset.BuildMomDiscOperator(options.scattering_order,
                                  options.geometry_type);
    groupset.BuildSubsets();

    ComputeSweepOrderings(groupset);
    InitFluxDataStructures(groupset);

    InitWGDSA(groupset);
    InitTGDSA(groupset);

    SolveGroupset(groupset, gs);

    CleanUpWGDSA(groupset);
    CleanUpTGDSA(groupset);

    ResetSweepOrderings(groupset);

    MPI_Barrier(MPI_COMM_WORLD);
  }

  chi_log.Log(LOG_0) << "NPTransport solver execution completed\n";
}


//###################################################################
/**Solves a single groupset.*/
void LinearBoltzmann::Solver::SolveGroupset(LBSGroupset& groupset,
                                            int group_set_num)
{
  source_event_tag = chi_log.GetRepeatingEventTag("Set Source");

  //================================================== Setting up required
  //                                                   sweep chunks
  SweepChunk* sweep_chunk = SetSweepChunk(groupset);
  MainSweepScheduler sweepScheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
                                    &groupset.angle_agg);

  *(char *)0 = 0;

  if (groupset.iterative_method == NPT_CLASSICRICHARDSON)
  {
    ClassicRichardson(groupset, group_set_num, sweep_chunk, sweepScheduler);
  }
  else if (groupset.iterative_method == NPT_GMRES)
  {
    GMRES(groupset, group_set_num, sweep_chunk, sweepScheduler);
  }

  delete sweep_chunk;

  if (options.write_restart_data)
    WriteRestartData(options.write_restart_folder_name,
                     options.write_restart_file_base);

  chi_log.Log(LOG_0)
    << "Groupset solve complete.                  Process memory = "
    << std::setprecision(3)
    << chi_console.GetMemoryUsageInMB() << " MB";
}

