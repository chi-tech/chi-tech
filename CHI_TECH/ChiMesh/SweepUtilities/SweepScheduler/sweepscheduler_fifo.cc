#include "sweepscheduler.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;


//###################################################################
/**Applies a First-In-First-Out sweep scheduling.*/
void chi_mesh::sweep_management::SweepScheduler::ScheduleAlgoFIFO()
{
  chi_log.LogEvent(sweep_event_tag, ChiLog::EventType::EVENT_BEGIN);

  auto ev_info_i =
    std::make_shared<ChiLog::EventInfo>(std::string("Sweep initiated"));

  chi_log.LogEvent(sweep_event_tag,
                   ChiLog::EventType::SINGLE_OCCURRENCE,ev_info_i);

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
        AngleSetGroupAdvance(sweep_chunk, q, sweep_timing_events_tag);
    }
  }

  //================================================== Reset all
  for (auto angsetgroup : angle_agg->angle_set_groups)
    angsetgroup->ResetSweep();

  for (auto bndry : angle_agg->sim_boundaries)
  {
    if (bndry->Type() == chi_mesh::sweep_management::BoundaryType::REFLECTING)
    {
      auto rbndry = (chi_mesh::sweep_management::BoundaryReflecting*)bndry;
      rbndry->ResetAnglesReadyStatus();
    }
  }

  //================================================== Receive delayed data
  MPI_Barrier(MPI_COMM_WORLD);
  for (auto sorted_angleset : rule_values)
  {
    TAngleSet *angleset = sorted_angleset.angle_set;
    angleset->ReceiveDelayedData(sorted_angleset.set_index);
  }

  chi_log.LogEvent(sweep_event_tag, ChiLog::EventType::EVENT_END);

}