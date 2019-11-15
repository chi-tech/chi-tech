#include "sweepscheduler.h"

//###################################################################
/**This is the entry point for sweeping.*/
void chi_mesh::sweep_management::SweepScheduler::
     Sweep(SweepChunk* in_sweep_chunk)
{
  sweep_chunk = in_sweep_chunk;

  if (scheduler_type == SchedulingAlgorithm::FIRST_IN_FIRST_OUT)
    ScheduleAlgoFIFO();
  else if (scheduler_type == SchedulingAlgorithm::DEPTH_OF_GRAPH)
    ScheduleAlgoDOG();
}