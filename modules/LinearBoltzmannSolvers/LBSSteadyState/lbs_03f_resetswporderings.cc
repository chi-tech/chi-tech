#include "lbs_linear_boltzmann_solver.h"

#include "ChiConsole/chi_console.h"

#include "ChiLog/chi_log.h"
;

#include "chi_mpi.h"
#include "LBSSteadyState/Groupset/lbs_groupset.h"


#include <iomanip>

//###################################################################
/**Clears all the sweep orderings for a groupset in preperation for
 * another.*/
void lbs::SteadyStateSolver::ResetSweepOrderings(LBSGroupset& groupset)
{
  chi::log.Log0Verbose1()
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

  if (options.verbose_inner_iterations)
    chi::log.Log()
      << "SPDS and FLUDS reset complete.            Process memory = "
      << std::setprecision(3)
      << chi_objects::ChiConsole::GetMemoryUsageInMB() << " MB";

  double local_app_memory =
    chi::log.ProcessEvent(chi_objects::ChiLog::StdTags::MAX_MEMORY_USAGE,
                          chi_objects::ChiLog::EventOperation::MAX_VALUE);
  double total_app_memory=0.0;
  MPI_Allreduce(&local_app_memory,&total_app_memory,
                1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  double max_proc_memory=0.0;
  MPI_Allreduce(&local_app_memory,&max_proc_memory,
                1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (options.verbose_inner_iterations)
    chi::log.Log()
      << "\n" << std::setprecision(3)
      << "           Total application memory (max): "
      << total_app_memory/1024.0 << " GB\n"
      << "           Maximum process memory        : "
      << max_proc_memory/1024.0 << " GB\n\n";

}