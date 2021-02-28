#include "lbs_linear_boltzmann_solver.h"

#include <ChiConsole/chi_console.h>
extern ChiConsole&  chi_console;

#include <chi_log.h>
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include <iomanip>

//###################################################################
/**Clears all the sweep orderings for a groupset in preperation for
 * another.*/
void LinearBoltzmann::Solver::ResetSweepOrderings(LBSGroupset& groupset)
{
  chi_log.Log(LOG_0VERBOSE_1)
    << "Resetting SPDS and FLUDS";

  groupset.sweep_orderings.clear();

  auto& angle_agg = groupset.angle_agg;

  for (auto& angset_grp : angle_agg.angle_set_groups)
  {
    for (auto& angset : angset_grp.angle_sets)
    {
      delete angset->fluds;
    }
    angset_grp.angle_sets.clear();
  }
  angle_agg.angle_set_groups.clear();

  MPI_Barrier(MPI_COMM_WORLD);

  chi_log.Log(LOG_0)
    << "SPDS and FLUDS reset complete.            Process memory = "
    << std::setprecision(3)
    << chi_console.GetMemoryUsageInMB() << " MB";

  double local_app_memory =
    chi_log.ProcessEvent(ChiLog::StdTags::MAX_MEMORY_USAGE,
                         ChiLog::EventOperation::MAX_VALUE);
  double total_app_memory=0.0;
  MPI_Allreduce(&local_app_memory,&total_app_memory,
                1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  double max_proc_memory=0.0;
  MPI_Allreduce(&local_app_memory,&max_proc_memory,
                1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  chi_log.Log(LOG_0)
    << "\n" << std::setprecision(3)
    << "           Total application memory (max): "
    << total_app_memory/1000.0 << " GB\n"
    << "           Maximum process memory        : "
    << max_proc_memory/1000.0 << " GB\n\n";

}