#include "sweepscheduler.h"

#include <chi_mpi.h>
#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

extern double chi_global_timings[20];

//###################################################################
/**Applies a First-In-First-Out sweep scheduling.*/
void chi_mesh::sweep_management::SweepScheduler::ScheduleAlgoFIFO()
{
  ChiTimer t16_sweeptime; t16_sweeptime.Reset();
  //================================================== Loop over AngleSetGroups
  // For 3D geometry this will be 8, one for each octant.
  // For 2D geometry this will be 4, one for each quadrant.
  // For 1D geometry this will be 2, one for left and one for right
  AngleSetStatus completion_status = AngleSetStatus::NOT_FINISHED;
  while (completion_status == AngleSetStatus::NOT_FINISHED)
  {
    completion_status = AngleSetStatus::FINISHED;
    for (int q=0; q<angle_agg->angle_set_groups.size(); q++)
    {
      completion_status = angle_agg->angle_set_groups[q]->
        AngleSetGroupAdvance(sweep_chunk, q);
    }
  }

  //================================================== Reset all
  for (auto angsetgroup : angle_agg->angle_set_groups)
    angsetgroup->ResetSweep();

  //================================================== Receive delayed data
  MPI_Barrier(MPI_COMM_WORLD);
  for (auto sorted_angleset : rule_values)
  {
    TAngleSet *angleset = sorted_angleset.angle_set;
    angleset->ReceiveDelayedData(sorted_angleset.set_index);
  }

  chi_global_timings[16] += t16_sweeptime.GetTime()/1000.0;
  chi_global_timings[17] += 1.0;

}