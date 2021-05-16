#include "lbs_linear_boltzmann_solver.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Sweeps all groupsets with the latest available phi_old.*/
void LinearBoltzmann::Solver::SweepAllGroupSets()
{
  MPI_Barrier(MPI_COMM_WORLD);
  int gs=-1;
  for (auto& groupset : group_sets)
  {
    ++gs;
    chi_log.Log(LOG_0)
      << "\n********* Initializing Groupset " << gs << "\n" << std::endl;

    ComputeSweepOrderings(groupset);
    InitFluxDataStructures(groupset);

    SweepChunk* sweep_chunk = SetSweepChunk(groupset);
    MainSweepScheduler sweepScheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
                                      &groupset.angle_agg);

    SetSource(groupset,SourceFlags::USE_MATERIAL_SOURCE,false);

    groupset.ZeroPsiDataStructures();
    phi_new_local.assign(phi_new_local.size(),0.0); //Ensure phi_new=0.0
    sweepScheduler.Sweep(sweep_chunk);

    delete sweep_chunk;

    ResetSweepOrderings(groupset);

    MPI_Barrier(MPI_COMM_WORLD);
  }

  chi_log.Log(LOG_0) << "All groupset sweep-execution completed\n";
}