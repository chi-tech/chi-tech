#include "chi_sweepscheduler.h"

#include <chi_mpi.h>
#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

extern CHI_MPI chi_mpi;
extern CHI_LOG chi_log;

extern double chi_global_timings[20];

//###################################################################
/**Applies a First-In-First-Out sweep scheduling.*/
void chi_mesh::SweepManagement::SweepScheduler::ScheduleAlgoFIFO()
{
  ChiTimer t16_sweeptime; t16_sweeptime.Reset();
  //================================================== Loop over AngleSetGroups
  // For 3D geometry this will be 8, one for each octant.
  // For 2D geometry this will be 4, one for each quadrant.
  // For 1D geometry this will be 2, one for left and one for right
  bool completion_status = FLAG_FINISHED;
  bool first_group = true;
  while ((completion_status == FLAG_NOT_FINISHED) || (first_group))
  {
    completion_status = FLAG_FINISHED;
    first_group=false;
    for (int q=0; q<angle_agg->angle_set_groups.size(); q++)
    {
      bool status = angle_agg->angle_set_groups[q]->
        AngleSetGroupAdvance(sweep_chunk, q);
      completion_status = completion_status && status;
    }
  }

  //================================================== Reset all
  for (int q=0; q<angle_agg->angle_set_groups.size(); q++)
  {
    angle_agg->angle_set_groups[q]->ResetSweep();
  }
  chi_global_timings[16] += t16_sweeptime.GetTime()/1000.0;
  chi_global_timings[17] += 1.0;

}