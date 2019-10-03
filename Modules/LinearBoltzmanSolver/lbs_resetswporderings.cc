#include "lbs_linear_boltzman_solver.h"

#include <ChiConsole/chi_console.h>
#include <chi_mpi.h>

extern ChiConsole chi_console;
extern ChiMPI chi_mpi;

extern double chi_global_timings[20];

//###################################################################
/**Clears all the sweep orderings for a groupset in preperation for
 * another.*/
void LinearBoltzman::Solver::ResetSweepOrderings(LBSGroupset *groupset)
{
  chi_log.Log(LOG_0VERBOSE_1)
    << "Resetting SPDS and FLUDS";
  for (int so=0; so<sweep_orderings.size(); so++)
  {
    chi_mesh::SweepManagement::SPDS* cur_so =
      sweep_orderings[so];

    //delete cur_so->spls->fluds;
    delete cur_so->spls;
    delete cur_so;
  }

  sweep_orderings.clear();

  chi_mesh::SweepManagement::AngleAggregation* angle_agg = groupset->angle_agg;

  for (int asg=0; asg<angle_agg->angle_set_groups.size(); asg++)
  {
    chi_mesh::SweepManagement::AngleSetGroup* angset_grp =
      angle_agg->angle_set_groups[asg];

    for (int as=0; as<angset_grp->angle_sets.size(); as++)
    {
      chi_mesh::SweepManagement::AngleSet* angset =
        angset_grp->angle_sets[as];

      angset->local_psi.clear();
      angset->local_psi.shrink_to_fit();

      angset->deplocI_outgoing_psi.clear();
      angset->deplocI_outgoing_psi.shrink_to_fit();

      angset->prelocI_outgoing_psi.clear();
      angset->prelocI_outgoing_psi.shrink_to_fit();


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

  double total_app_memory=0.0;
  MPI_Allreduce(&chi_global_timings[9],&total_app_memory,
                1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  double max_proc_memory=0.0;
  MPI_Allreduce(&chi_global_timings[9],&max_proc_memory,
                1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  double local_avg_memory= chi_global_timings[10]/chi_global_timings[11];
  double avg_proc_memory=0.0;
  MPI_Allreduce(&local_avg_memory,&avg_proc_memory,
                1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  avg_proc_memory/=chi_mpi.process_count;
  chi_log.Log(LOG_0)
    << "\n" << std::setprecision(3)
    << "           Total application memory (max): " << total_app_memory/1000.0 << " GB\n"
    << "           Maximum process memory        : " << max_proc_memory/1000.0 << " GB\n"
    << "           Average process memory        : " << avg_proc_memory/1000.0 << " GB\n\n";

}