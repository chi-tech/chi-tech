#include "sweepscheduler.h"

#include <chi_log.h>
extern ChiLog chi_log;

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

//###################################################################
/**Get average sweep time from logging system.*/
double chi_mesh::sweep_management::SweepScheduler::GetAverageSweepTime()
{
  return chi_log.ProcessEvent(sweep_event_tag,
                              ChiLog::EventOperation::AVERAGE_DURATION);
}