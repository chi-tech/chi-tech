#include "chi_sweepscheduler.h"

//###################################################################
/**This is the entry point for sweeping.*/
void chi_mesh::SweepManagement::SweepScheduler::
     Sweep(SweepChunk* in_sweep_chunk)
{
  sweep_chunk = in_sweep_chunk;

  if (scheduler_type == FIFO)
    ScheduleAlgoFIFO();
  else if (scheduler_type == DEPTH_OF_GRAPH)
    ScheduleAlgoDOG();
}