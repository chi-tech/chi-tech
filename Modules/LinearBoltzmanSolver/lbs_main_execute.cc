#include "lbs_linear_boltzman_solver.h"
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
void LinearBoltzman::Solver::Execute()
{
  MPI_Barrier(MPI_COMM_WORLD);
  for (int gs=0; gs<group_sets.size(); gs++)
  {
    chi_log.Log(LOG_0)
      << "\n********* Initializing Groupset " << gs << "\n" << std::endl;

    group_sets[gs]->BuildDiscMomOperator(options.scattering_order,
                                         options.geometry_type);
    group_sets[gs]->BuildMomDiscOperator(options.scattering_order,
                                         options.geometry_type);
    group_sets[gs]->BuildSubsets();

    ComputeSweepOrderings(group_sets[gs]);
    InitFluxDataStructures(group_sets[gs]);

    InitWGDSA(group_sets[gs]);
    InitTGDSA(group_sets[gs]);

    SolveGroupset(gs);

    CleanUpWGDSA(group_sets[gs]);
    CleanUpTGDSA(group_sets[gs]);

    ResetSweepOrderings(group_sets[gs]);

    MPI_Barrier(MPI_COMM_WORLD);
  }



  chi_log.Log(LOG_0) << "NPTransport solver execution completed\n";
}


//###################################################################
/**Solves a single groupset.*/
void LinearBoltzman::Solver::SolveGroupset(int group_set_num)
{
  source_event_tag = chi_log.GetRepeatingEventTag("Set Source");
  LBSGroupset* group_set = group_sets[group_set_num];

  //================================================== Setting up required
  //                                                   sweep chunks
  SweepChunk* sweep_chunk = SetSweepChunk(group_set_num);
  MainSweepScheduler sweepScheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
                                    group_set->angle_agg);

  if (group_set->iterative_method == NPT_CLASSICRICHARDSON)
  {
    ClassicRichardson(group_set_num, sweep_chunk, sweepScheduler);
  }
  else if (group_set->iterative_method == NPT_GMRES)
  {
    GMRES(group_set_num, sweep_chunk, sweepScheduler);
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

