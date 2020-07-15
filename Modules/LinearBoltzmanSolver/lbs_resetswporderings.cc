#include "lbs_linear_boltzman_solver.h"

#include <ChiConsole/chi_console.h>
#include <chi_mpi.h>

extern ChiConsole&  chi_console;
extern ChiMPI& chi_mpi;


#include <chi_log.h>

extern ChiLog& chi_log;

#include <iomanip>

//###################################################################
/**Clears all the sweep orderings for a groupset in preperation for
 * another.*/
void LinearBoltzman::Solver::ResetSweepOrderings(LBSGroupset *groupset)
{
  chi_log.Log(LOG_0VERBOSE_1)
    << "Resetting SPDS and FLUDS";
  for (int so=0; so<sweep_orderings.size(); so++)
  {
    chi_mesh::sweep_management::SPDS* cur_so =
      sweep_orderings[so];

    //delete cur_so->spls->fluds;
    delete cur_so->spls;
    delete cur_so;
  }

  sweep_orderings.clear();

  chi_mesh::sweep_management::AngleAggregation* angle_agg = groupset->angle_agg;

  for (int asg=0; asg<angle_agg->angle_set_groups.size(); asg++)
  {
    chi_mesh::sweep_management::AngleSetGroup* angset_grp =
      angle_agg->angle_set_groups[asg];

    for (int as=0; as<angset_grp->angle_sets.size(); as++)
    {
      chi_mesh::sweep_management::AngleSet* angset =
        angset_grp->angle_sets[as];

      delete angset->fluds;
      delete angset;
    }
    angset_grp->angle_sets.clear();
    delete angset_grp;
  }
  angle_agg->angle_set_groups.clear();
  delete angle_agg;

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