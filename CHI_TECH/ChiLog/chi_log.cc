#include "chi_log.h"
#include <chi_mpi.h>
#include <ChiTimer/chi_timer.h>

extern ChiMPI    chi_mpi;
extern ChiTimer  chi_program_timer;

#include <sstream>
#include <iomanip>

//###################################################################
/** Default constructor*/
ChiLog::ChiLog() noexcept
{
  verbosity = LOG_0VERBOSE_0;
  std::string memory_usage_event("Maximum Memory Usage");
  repeating_events.emplace_back(memory_usage_event);

  RepeatingEvent& ref_rep_event = repeating_events.back();

  ref_rep_event.events.emplace_back(
    chi_program_timer.GetTime(),
    EventType::EVENT_CREATED,
    std::make_shared<EventInfo>());
}

//###################################################################
/** Makes a log entry.*/
LogStream ChiLog::Log(LOG_LVL level)
{
  switch (level)
  {
    case LOG_0:
    {
      if (chi_mpi.location_id == 0)
      {
        std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
        return LogStream(&std::cout, header);
      }
      else
      {
        std::string header = " ";
        return LogStream(&dummy_stream, header);
      }
    }
    case LOG_0WARNING:
    {
      if (chi_mpi.location_id == 0)
      {
        std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
        header += "**WARNING** ";
        return LogStream(&std::cout, header);
      }
      else
      {
        std::string header = " ";
        return LogStream(&dummy_stream, header);
      }
    }
    case LOG_0ERROR:
    {
      if (chi_mpi.location_id == 0)
      {
        std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
        header += "**!**ERROR**!** ";
        return LogStream(&std::cerr, header);
      }
      else
      {
        std::string header = " ";
        return LogStream(&dummy_stream, header);
      }
    }
    case LOG_0VERBOSE_0:
    case LOG_0VERBOSE_1:
    case LOG_0VERBOSE_2:
    {
      if ((chi_mpi.location_id == 0) && (verbosity >= level))
      {
        std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
        return LogStream(&std::cout, header);
      }
      else
      {
        std::string header = " ";
        return LogStream(&dummy_stream, header);
      }
    }
    case LOG_ALL:
    {
      std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
      return LogStream(&std::cout, header);
    }
    case LOG_ALLWARNING:
    {
      std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
      header += "**WARNING** ";
      return LogStream(&std::cout, header);
    }
    case LOG_ALLERROR:
    {
      std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
      header += "**!**ERROR**!** ";
      return LogStream(&std::cerr, header);
    }

    case LOG_ALLVERBOSE_0:
    case LOG_ALLVERBOSE_1:
    case LOG_ALLVERBOSE_2:
    {
      if (verbosity >= (level-6))
      {
        std::string header = "[" + std::to_string(chi_mpi.location_id) + "]  ";
        return LogStream(&std::cout, header);
      }
      else
      {
        std::string header = " ";
        return LogStream(&dummy_stream, header);
      }
    }
  }


  std::string header = " ";
  return LogStream(&dummy_stream, header);
}


//###################################################################
/** Sets the verbosity level.*/
void ChiLog::SetVerbosity(int int_level)
{
  if (int_level == 0)
  {
    verbosity = LOG_0VERBOSE_0;
  }
  else if (int_level == 1)
  {
    verbosity = LOG_0VERBOSE_1;
  }
  else if (int_level == 2)
  {
    verbosity = LOG_0VERBOSE_2;
  }
}

//###################################################################
/** Gets the current verbosity level.*/
int ChiLog::GetVerbosity()
{
  return verbosity;
}

//###################################################################
/** Returns a unique tag to a newly created repeating event.*/
size_t ChiLog::GetRepeatingEventTag(std::string event_name)
{
  repeating_events.emplace_back(event_name);

  RepeatingEvent& ref_rep_event = repeating_events.back();

  ref_rep_event.events.emplace_back(
    chi_program_timer.GetTime(),
    EventType::EVENT_CREATED,
    std::make_shared<EventInfo>());

  return repeating_events.size()-1;
}

//###################################################################
/**Logs an event with the supplied event information.*/
void ChiLog::LogEvent(size_t ev_tag,
                      ChiLog::EventType ev_type,
                      std::shared_ptr<ChiLog::EventInfo> ev_info)
{
  if (ev_tag >= repeating_events.size())
    return;

  RepeatingEvent& ref_rep_event = repeating_events[ev_tag];

  ref_rep_event.events.emplace_back(
    chi_program_timer.GetTime(),
    ev_type,
    ev_info);
}

//###################################################################
/**Logs an event without any event information.*/
void ChiLog::LogEvent(size_t ev_tag,
                      ChiLog::EventType ev_type)
{
  if (ev_tag >= repeating_events.size())
    return;

  RepeatingEvent& ref_rep_event = repeating_events[ev_tag];

  ref_rep_event.events.emplace_back(
    chi_program_timer.GetTime(),
    ev_type,
    nullptr);
}

//###################################################################
/**Returns a string representation of the event history associated with
 * the tag. Each event entry will be prepended by the location id and
 * the program timestamp in seconds. This method uses the
 * ChiLog::EventInfo::GetString method to append information. This allows
 * derived classes to implement more sophisticated outputs.*/
std::string ChiLog::PrintEventHistory(size_t ev_tag)
{
  std::stringstream outstr;
  if (ev_tag >= repeating_events.size())
    return outstr.str();

  RepeatingEvent& ref_rep_event = repeating_events[ev_tag];

  for (auto& event : ref_rep_event.events)
  {
    outstr << "[" << chi_mpi.location_id << "] ";

    char buf[100];
    sprintf(buf,"%16.9f",event.ev_time/1000.0);
    outstr << buf << " ";

    switch (event.ev_type)
    {
      case EventType::EVENT_CREATED:
        outstr << "EVENT_CREATED ";
        break;
      case EventType::SINGLE_OCCURRENCE:
        outstr << "SINGLE_OCCURRENCE ";
        break;
      case EventType::EVENT_BEGIN:
        outstr << "EVENT_BEGIN ";
        break;
      case EventType::EVENT_END:
        outstr << "EVENT_END ";
        break;
    }

    if (event.ev_info != nullptr)
      outstr << event.ev_info->GetString();
    outstr << std::endl;
  }

  return outstr.str();
}

//###################################################################
/**Processes an event given an event operation. See ChiLog for further
 * reference.*/
double ChiLog::ProcessEvent(size_t ev_tag,
                            ChiLog::EventOperation ev_operation)
{
  if (ev_tag >= repeating_events.size())
    return 0.0;

  RepeatingEvent& ref_rep_event = repeating_events[ev_tag];

  double ret_val = 0.0;
  switch (ev_operation)
  {
    case EventOperation::NUMBER_OF_OCCURRENCES:
    {
      for (auto& event : ref_rep_event.events)
      {
        if ( (event.ev_type == EventType::EVENT_CREATED) or
             (event.ev_type == EventType::SINGLE_OCCURRENCE) or
             (event.ev_type == EventType::EVENT_BEGIN))
          ret_val += 1.0;
      }//for events
      break;
    }
    case EventOperation::TOTAL_DURATION:
    {
      double start_time = 0.0;
      for (auto& event : ref_rep_event.events)
      {
        if (event.ev_type == EventType::EVENT_BEGIN)
          start_time = event.ev_time;
        if (event.ev_type == EventType::EVENT_END)
          ret_val += event.ev_time - start_time;
      }//for events
      ret_val *= 1000.0;
      break;
    }
    case EventOperation::AVERAGE_DURATION:
    {
      double start_time = 0.0;
      int counter = 0;
      for (auto& event : ref_rep_event.events)
      {
        if (event.ev_type == EventType::EVENT_BEGIN)
          start_time = event.ev_time;

        if (event.ev_type == EventType::EVENT_END)
        {
          ret_val += event.ev_time - start_time;
          counter++;
        }
      }//for events
      ret_val /= (1000.0*counter);
      break;
    }
    case EventOperation::MAX_VALUE:
    {
      ret_val = 0.0;
      for (auto& event : ref_rep_event.events)
      {
        if ((event.ev_type == EventType::SINGLE_OCCURRENCE) or
            (event.ev_type == EventType::EVENT_BEGIN) or
            (event.ev_type == EventType::EVENT_END))
        {
          if (event.ev_info != nullptr)
            ret_val = std::max(event.ev_info->arb_value,ret_val);
        }
      }//for events
      break;
    }
    case EventOperation::AVERAGE_VALUE:
    {
      ret_val = 0.0;
      int count = 0;
      for (auto& event : ref_rep_event.events)
      {
        if ((event.ev_type == EventType::SINGLE_OCCURRENCE) or
            (event.ev_type == EventType::EVENT_BEGIN) or
            (event.ev_type == EventType::EVENT_END))
        {
          if (event.ev_info != nullptr)
          {
            ret_val += event.ev_info->arb_value;
            ++count;
          }
        }
      }//for events
      if (count == 0) count = 1;
      ret_val /= count;
      break;
    }
  }//switch

  return ret_val;
}
