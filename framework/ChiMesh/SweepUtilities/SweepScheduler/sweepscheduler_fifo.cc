#include "sweepscheduler.h"

#include "chi_mpi.h"
#include "chi_log.h"

//###################################################################
/**Applies a First-In-First-Out sweep scheduling.*/
void chi_mesh::sweep_management::SweepScheduler::
  ScheduleAlgoFIFO(SweepChunk& sweep_chunk)
{
  typedef AngleSetStatus Status;

  chi::log.LogEvent(sweep_event_tag, chi_objects::ChiLog::EventType::EVENT_BEGIN);

  auto ev_info_i =
    std::make_shared<chi_objects::ChiLog::EventInfo>(std::string("Sweep initiated"));

  chi::log.LogEvent(sweep_event_tag,
                   chi_objects::ChiLog::EventType::SINGLE_OCCURRENCE, ev_info_i);

  //================================================== Loop over AngleSetGroups
  AngleSetStatus completion_status = AngleSetStatus::NOT_FINISHED;
  while (completion_status == AngleSetStatus::NOT_FINISHED)
  {
    completion_status = AngleSetStatus::FINISHED;
    for (int q=0; q<angle_agg.angle_set_groups.size(); q++)
    {
      completion_status = angle_agg.angle_set_groups[q].
        AngleSetGroupAdvance(sweep_chunk, q, sweep_timing_events_tag);
    }
  }

  //================================================== Receive delayed data
  MPI_Barrier(MPI_COMM_WORLD);
  bool received_delayed_data = false;
  while (not received_delayed_data)
  {
    received_delayed_data = true;
    for (auto& sorted_angleset : rule_values)
    {
      auto& as = sorted_angleset.angle_set;

      if (as->FlushSendBuffers() == Status::MESSAGES_PENDING)
        received_delayed_data = false;

      if (not as->ReceiveDelayedData(sorted_angleset.set_index))
        received_delayed_data = false;
    }
  }

  //================================================== Reset all
  for (auto& angsetgroup : angle_agg.angle_set_groups)
    angsetgroup.ResetSweep();

  for (auto& [bid, bndry] : angle_agg.sim_boundaries)
  {
    if (bndry->Type() == chi_mesh::sweep_management::BoundaryType::REFLECTING)
    {
      auto rbndry = std::static_pointer_cast<
        chi_mesh::sweep_management::BoundaryReflecting>(bndry);
      rbndry->ResetAnglesReadyStatus();
    }
  }

//  //================================================== Receive delayed data
//  MPI_Barrier(MPI_COMM_WORLD);
//  for (auto& sorted_angleset : rule_values)
//  {
//    auto angleset = sorted_angleset.angle_set;
//    angleset->ReceiveDelayedData(sorted_angleset.set_index);
//  }

  chi::log.LogEvent(sweep_event_tag, chi_objects::ChiLog::EventType::EVENT_END);

}