#include "sweepscheduler.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Sweep scheduler constructor*/
chi_mesh::sweep_management::SweepScheduler::SweepScheduler(
    SchedulingAlgorithm in_scheduler_type,
    chi_mesh::sweep_management::AngleAggregation& in_angle_agg,
    SweepChunk& in_sweep_chunk) :
  scheduler_type(in_scheduler_type),
  angle_agg(in_angle_agg),
  m_sweep_chunk(in_sweep_chunk),
  sweep_timing_events_tag({Chi::log.GetRepeatingEventTag("Sweep Chunk Only Timing")
  }),
  sweep_event_tag(Chi::log.GetRepeatingEventTag("Sweep Timing"))
{
  angle_agg.InitializeReflectingBCs();

  if (scheduler_type == SchedulingAlgorithm::DEPTH_OF_GRAPH)
    InitializeAlgoDOG();

  //=================================== Initialize delayed upstream data
  for (auto& angsetgrp : in_angle_agg.angle_set_groups)
    for (auto& angset : angsetgrp.angle_sets)
      angset->InitializeDelayedUpstreamData();

  //=================================== Get local max num messages accross
  //                                    anglesets
  int local_max_num_messages = 0;
  for (auto& angsetgrp : in_angle_agg.angle_set_groups)
    for (auto& angset : angsetgrp.angle_sets)
      local_max_num_messages = std::max(angset->GetMaxBufferMessages(),
                                        local_max_num_messages);

  //=================================== Reconcile all local maximums
  int global_max_num_messages = 0;
  MPI_Allreduce(&local_max_num_messages,
                &global_max_num_messages,
                1, MPI_INT,
                MPI_MAX, Chi::mpi.comm);

  //=================================== Propogate items back to sweep buffers
  for (auto& angsetgrp : in_angle_agg.angle_set_groups)
    for (auto& angset : angsetgrp.angle_sets)
      angset->SetMaxBufferMessages(global_max_num_messages);
}

//###################################################################
/**Returns the referenced sweep chunk.*/
chi_mesh::sweep_management::SweepChunk&
  chi_mesh::sweep_management::SweepScheduler::GetSweepChunk()
{
  return m_sweep_chunk;
}