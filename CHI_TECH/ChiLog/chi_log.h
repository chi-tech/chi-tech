#ifndef _chi_log_h
#define _chi_log_h

#include "chi_logstream.h"
#include <vector>
#include <memory>

/**Logging level*/
enum LOG_LVL {LOG_0=1,                //Used only for location 0
              LOG_0WARNING=2,         //Warning only for location 0
              LOG_0ERROR=3,           //Error only for location 0
              LOG_0VERBOSE_0=4,       //Default verbosity level
              LOG_0VERBOSE_1=5,       //Used only if verbosity level equals 1
              LOG_0VERBOSE_2=6,       //Used only if verbosity level equals 2
              LOG_ALL=7,              //Verbose level 0 all locations
              LOG_ALLWARNING=8,       //Warning for any location
              LOG_ALLERROR=9,         //Error for any location
              LOG_ALLVERBOSE_0=10,    //Default verbosity level
              LOG_ALLVERBOSE_1=11,    //Used only if verbosity level equals 1
              LOG_ALLVERBOSE_2=12};   //Used only if verbosity level equals 2



//###################################################################
/**Object for controlling logging.
 *
 * Part A: Output logs
 * There are three levels of
 * verbosity in ChiTech: Zero(Default), One and Two. These can
 * be set on the command line via the switch -v followed by a
 * space and the number for the verbosity (0,1 or 2). The lua command
 * chiLogSetVerbosity(int_level) achieves the same.
 *
 * Part B: Repeating events log
 * Suppose a routine gets called multiple times. Now suppose we want to know
 * how many times this routine was called, how long it ran (total and on-avg),
 * and we want to know at which timestamp it ran every time. We do this
 * by using a repeating event.*/
class ChiLog
{
private:
  DummyStream dummy_stream;
  int verbosity;

public:
  //00
                  ChiLog() noexcept;
  //01
  LogStream       Log(LOG_LVL level);
  void            SetVerbosity(int int_level);
  int             GetVerbosity();

public:
  class RepeatingEvent;
  enum StdTags
  {
    MAX_MEMORY_USAGE = 0
  };
  enum class EventType
  {
    EVENT_CREATED = 0,
    SINGLE_OCCURRENCE = 1,
    EVENT_BEGIN = 2,
    EVENT_END = 3
  };
  enum class EventOperation
  {
    NUMBER_OF_OCCURRENCES = 0,
    TOTAL_DURATION = 1,
    AVERAGE_DURATION = 2,
    MAX_VALUE = 3
  };
  struct EventInfo;
  struct Event;

private:
  std::vector<RepeatingEvent> repeating_events;

public:
  size_t GetRepeatingEventTag(std::string event_name);
  void   LogEvent(size_t ev_tag,
                  EventType ev_type,
                  std::shared_ptr<EventInfo>& ev_info);
  void   LogEvent(size_t ev_tag,
                  EventType ev_type);
  void   PrintEventHistory(size_t ev_tag, std::ostream& outstr);
  double ProcessEvent(size_t ev_tag, EventOperation ev_operation);
};

//###################################################################
/** */
struct ChiLog::EventInfo
{
  double arb_value = 0.0;
           EventInfo()= default;
  explicit EventInfo(double in_value) : arb_value(in_value) {}

  virtual std::string GetString()
  {
    return std::to_string(arb_value);
  }
};

//###################################################################
/** */
struct ChiLog::Event
{
  const double     ev_time = 0.0;
  const EventType  ev_type = EventType::SINGLE_OCCURRENCE;
  std::shared_ptr<EventInfo> ev_info;

  Event(double in_time,
        EventType in_ev_type,
        std::shared_ptr<EventInfo> in_event_info) :
    ev_time(in_time),
    ev_type(in_ev_type),
    ev_info(in_event_info)
  {}
};

//###################################################################
/**Repeating event object.*/
class ChiLog::RepeatingEvent
{
public:
  const std::string  name;
  std::vector<Event> events;
public:
  explicit RepeatingEvent(std::string& in_event_name) : name(in_event_name)
  {  }
};

#endif