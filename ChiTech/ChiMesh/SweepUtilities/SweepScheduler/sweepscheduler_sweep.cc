#include "sweepscheduler.h"

#include <chi_log.h>
extern ChiLog& chi_log;

//###################################################################
/**This is the entry point for sweeping.*/
void chi_mesh::sweep_management::SweepScheduler::
     Sweep(SweepChunk& sweep_chunk)
{
  if (scheduler_type == SchedulingAlgorithm::FIRST_IN_FIRST_OUT)
    ScheduleAlgoFIFO(sweep_chunk);
  else if (scheduler_type == SchedulingAlgorithm::DEPTH_OF_GRAPH)
    ScheduleAlgoDOG(sweep_chunk);
}

//###################################################################
/**Get average sweep time from logging system.*/
double chi_mesh::sweep_management::SweepScheduler::GetAverageSweepTime()
{
  return chi_log.ProcessEvent(sweep_event_tag,
                              ChiLog::EventOperation::AVERAGE_DURATION);
}

//###################################################################
/**Get relevant sweep timing information.
 *
 * [0] Total sweep time
 * [1] Total chunk time
 * [2] Total chunk time / total sweep time
 * */
std::vector<double>
  chi_mesh::sweep_management::SweepScheduler::GetAngleSetTimings()
{
  std::vector<double> info;

  double total_sweep_time =
    chi_log.ProcessEvent(sweep_event_tag,
                         ChiLog::EventOperation::TOTAL_DURATION);

  double total_chunk_time =
    chi_log.ProcessEvent(sweep_timing_events_tag[0],
                         ChiLog::EventOperation::TOTAL_DURATION);

  double ratio_sweep_to_chunk = total_chunk_time/total_sweep_time;

  info.push_back(total_sweep_time);
  info.push_back(total_chunk_time);
  info.push_back(ratio_sweep_to_chunk);

  return info;
}