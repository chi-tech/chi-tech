#include "sweepscheduler.h"

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Sweep scheduler constructor*/
chi_mesh::sweep_management::SweepScheduler::SweepScheduler(
    SchedulingAlgorithm in_scheduler_type,
    chi_mesh::sweep_management::AngleAggregation *in_angle_agg) :
  sweep_event_tag(chi_log.GetRepeatingEventTag("Sweep Timing")),
  sweep_timing_events_tag({
    chi_log.GetRepeatingEventTag("Sweep Chunk Only Timing")
  })
{
  scheduler_type = in_scheduler_type;
  angle_agg      = in_angle_agg;

  angle_agg->InitializeReflectingBCs();

  if (scheduler_type == SchedulingAlgorithm::DEPTH_OF_GRAPH)
    InitializeAlgoDOG();

  //=================================== Initialize delayed upstream data
  for (auto angsetgrp : in_angle_agg->angle_set_groups)
    for (auto angset : angsetgrp->angle_sets)
      angset->InitializeDelayedUpstreamData();

  //=================================== Get local max num messages accross
  //                                    anglesets
  int local_max_num_messages = 0;
  for (auto angsetgrp : in_angle_agg->angle_set_groups)
  {
    for (auto angset : angsetgrp->angle_sets)
    {
      local_max_num_messages = std::max(
        angset->GetMaxBufferMessages(),
        local_max_num_messages);
    }
  }

  //=================================== Reconcile all local maximums
  int global_max_num_messages = 0;
  MPI_Allreduce(&local_max_num_messages,
                &global_max_num_messages,
                1, MPI_INT,
                MPI_MAX, MPI_COMM_WORLD);

  //=================================== Propogate items back to sweep buffers
  for (auto angsetgrp : in_angle_agg->angle_set_groups)
    for (auto angset : angsetgrp->angle_sets)
      angset->SetMaxBufferMessages(global_max_num_messages);
}