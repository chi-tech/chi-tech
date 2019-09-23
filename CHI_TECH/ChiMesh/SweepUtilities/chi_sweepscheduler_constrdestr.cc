#include "chi_sweepscheduler.h"

//###################################################################
/**Sweep scheduler constructor*/
chi_mesh::SweepManagement::SweepScheduler::SweepScheduler(
    int in_scheduler_type,
    chi_mesh::SweepManagement::AngleAggregation *in_angle_agg)
{
  scheduler_type = in_scheduler_type;
  angle_agg      = in_angle_agg;

  if (scheduler_type == DEPTH_OF_GRAPH)
    InitializeAlgoDOG();

  for (auto angsetgrp : in_angle_agg->angle_set_groups)
  {
    for (auto angset : angsetgrp->angle_sets)
      angset->InitializeDelayedUpstreamData();
  }
}