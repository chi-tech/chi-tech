#ifndef _chi_log_h
#define _chi_log_h

#include "chi_logstream.h"
#include <vector>
#include <memory>

/**Logging level*/
enum LOG_LVL {LOG_0=1,                ///< Used only for location 0
              LOG_0WARNING=2,         ///< Warning only for location 0
              LOG_0ERROR=3,           ///< Error only for location 0
              LOG_0VERBOSE_0=4,       ///< Default verbosity level
              LOG_0VERBOSE_1=5,       ///< Used only if verbosity level equals 1
              LOG_0VERBOSE_2=6,       ///< Used only if verbosity level equals 2
              LOG_ALL=7,              ///< Verbose level 0 all locations
              LOG_ALLWARNING=8,       ///< Warning for any location
              LOG_ALLERROR=9,         ///< Error for any location
              LOG_ALLVERBOSE_0=10,    ///< Default verbosity level
              LOG_ALLVERBOSE_1=11,    ///< Used only if verbosity level equals 1
              LOG_ALLVERBOSE_2=12};   ///< Used only if verbosity level equals 2



//###################################################################
/**Object for controlling logging.
 *
 * ## Part A: Output logs
 * There are three levels of verbosity in ChiTech: Zero(Default), One and Two.
 * These can be set on the command line via the switch -v followed by a
 * space and the number for the verbosity (0,1 or 2).
 *
 * \code
 * ./bin/ChiTech InputFile.lua -v 1
 * \endcode
 *
 * The lua command `chiLogSetVerbosity(int_level)` achieves the same.\n\n
 *
 * Printing a log under the auspices of a verbosity level again has
 * numerous options. Firstly, any log can be a normal log, a warning
 * or an error. Secondly, a log can be either location 0 or all the locations
 * in a parallel process environment. The log option enums defined under
 * LOG_LVL are
 - LOG_0,                      Used only for location 0
 - LOG_0WARNING,               Warning only for location 0
 - LOG_0ERROR,                 Error only for location 0
 - LOG_0VERBOSE_0,             Default verbosity level
 - LOG_0VERBOSE_1,             Used only if verbosity level equals 1
 - LOG_0VERBOSE_2,             Used only if verbosity level equals 2
 - LOG_ALL,                    Verbose level 0 all locations
 - LOG_ALLWARNING,             Warning for any location
 - LOG_ALLERROR,               Error for any location
 - LOG_ALLVERBOSE_0,     Default verbosity level
 - LOG_ALLVERBOSE_1,     Used only if verbosity level equals 1
 - LOG_ALLVERBOSE_2,     Used only if verbosity level equals 2

 * A log can be made by first connecting code with the logger. This is done
 * by including the log header and then defining an extern reference to the
 * global object.
 *
 * \code
 * #include <chi_log.h>
 * extern ChiLog& chi_log;
 * \endcode
 *
 * A log is then written inside a piece of code as follows:
 *
 * \code
 * void PrintSomethingToLog()
 * {
 *   chi_log.Log(LOG_0) << "This is printed on location 0 only";
 *   chi_log.Log(LOG_0WARNING) << "This is a warning";
 *   chi_log.Log(LOG_0ERROR) << "This is an error";
 * }
 * \endcode
 *
\verbatim
[0]  This is printed on location 0 only
[0]  **WARNING** This is a warning
[0]  **!**ERROR**!** This is an error
\endverbatim
 *
 * ## Part B: Repeating events log
 * Suppose a routine or segment of code gets called multiple times. Now suppose
 * we want to know how many times this routine was called, how long it ran
 * (total and on-avg), or we want to know at which timestamp it ran every time.
 * We do this by using a repeating event. A repeating event is identified
 * by a tag. Therefore the first step is to obtain a unique id for the event using
 * ChiLog::GetRepeatingEventTag.
 *
 * \code
 * size_t tag = chi_log.GetRepeatingEventTag(std::string("Sweep Timing"))
 * \endcode
 *
 * Events can now be logged using
 * ChiLog::LogEvent(size_t ev_tag, ChiLog::EventType ev_type). As an example,
 * consider we want to count the number of occurrences of an event. We trigger
 * multiple occurences using
 *
 * \code
 * chi_log.LogEvent(tag,ChiLog::EventType::SINGLE_OCCURRENCE);
 * \endcode
 *
 * At the end of the block for which you wanted to log these events we can now
 * extract the number of occurrences using
 *
 * \code
 * double number_of_occ = chi_log.ProcessEvent(tag,ChiLog::EventOperation::NUMBER_OF_OCCURRENCES);
 *\endcode
 *
 * Below is a complete example:
 *
\code
size_t tag = chi_log.GetRepeatingEventTag(std::string());

for (int i=0; i<5; ++i)
  chi_log.LogEvent(tag,ChiLog::EventType::SINGLE_OCCURRENCE);

chi_log.Log(LOG_ALL)
  << chi_log.ProcessEvent(tag, ChiLog::EventOperation::NUMBER_OF_OCCURRENCES);
\endcode
\verbatim
6
\endverbatim

 * Since chi_log is a global object this functionality is really powerful for
 * use in modules or class objects where multiple methods can contribute to
 * the same event track. For example an object can have an event tag as a member
 * and initialize it at construction then all of the objects methods can
 * contribute to the event.
 *
 * ### Supplying event information
 * In addition to the ChiLog::EventType the user can also supply a reference to
 * a ChiLog::EventInfo structure. Be cautious with this though because each
 * event stored a double and a string which can look like a memory leak if
 * multiple events are pushed up with event information. Developers can supply
 * either a double or a string or both to an event info constructor to
 * instantiate an instance. The event arb_value is by default 0.0 and the event
 * arb_info is by default an empty string. An example is shown below:
 *
\code
size_t tag = chi_log.GetRepeatingEventTag(std::string());

for (int i=0; i<5; ++i)
{
  auto ev_info = std::make_shared<ChiLog::EventInfo>(i*2.0);
  chi_log.LogEvent(tag,ChiLog::EventType::SINGLE_OCCURRENCE,ev_info);
}

chi_log.Log(LOG_ALL)
  << chi_log.ProcessEvent(tag, ChiLog::EventOperation::AVERAGE_VALUE);
\endcode
\verbatim
4.0
\endverbatim
 *
 * Event information can also be supplied by string values but then the
 * average value is meaningless unless the average value is also supplied.
 * For example:
 *
\code
size_t tag = chi_log.GetRepeatingEventTag(std::string());

chi_log.LogEvent(tag,
                 ChiLog::EventType::SINGLE_OCCURRENCE,
                 std::make_shared<ChiLog::EventInfo>(std::string("A"),2.0));
chi_log.LogEvent(tag,
                 ChiLog::EventType::SINGLE_OCCURRENCE,
                 std::make_shared<ChiLog::EventInfo>(std::string("B")));
chi_log.LogEvent(tag,
                 ChiLog::EventType::SINGLE_OCCURRENCE,
                 std::make_shared<ChiLog::EventInfo>(std::string("C"),2.0));

chi_log.Log(LOG_ALL)
  << chi_log.ProcessEvent(tag, ChiLog::EventOperation::AVERAGE_VALUE);
\endcode
\verbatim
1.33333
\endverbatim
 *
 * To get a string value of the event history developers can use
 * ChiLog::PrintEventHistory along with the event tag. Just note that it will
 * automatically be formatted for each location so no need to use chi_log to
 * print it. Also, each event will be prepended with a program timestamp
 * in seconds.
 *
\code
std::cout << chi_log.PrintEventHistory(tag);
\endcode
\verbatim
1.33333
[0]      3.813119000 EVENT_CREATED
[0]      3.813120000 SINGLE_OCCURRENCE A
[0]      3.813121000 SINGLE_OCCURRENCE B
[0]      3.813122000 SINGLE_OCCURRENCE C
\endverbatim
 * */
class ChiLog
{
private:
  DummyStream dummy_stream;
  int verbosity;

  static ChiLog instance;

public:
  static ChiLog& GetInstance() noexcept
  { return instance;}
private:
  //00
                  ChiLog() noexcept;

public:
  //01
  LogStream       Log(LOG_LVL level=LOG_0);
  void            SetVerbosity(int int_level);
  int             GetVerbosity() const;

private:
  class RepeatingEvent;
public:
  enum StdTags
  {
    MAX_MEMORY_USAGE = 0 ///< Tag reserved for logging process memory
  };
  enum class EventType
  {
    EVENT_CREATED = 0,     ///< Automatically signalled when event is created
    SINGLE_OCCURRENCE = 1, ///< Signals a single occurrence
    EVENT_BEGIN = 2,       ///< Signals the begin of an event
    EVENT_END = 3          ///< Signals the end of an event
  };
  enum class EventOperation
  {
    NUMBER_OF_OCCURRENCES = 0, ///< Counts creation events, single occurrences and begins
    TOTAL_DURATION = 1,        ///< Integrates times between begins and ends
    AVERAGE_DURATION = 2,      ///< Computes average time between begins and ends
    MAX_VALUE = 3,             ///< Computes the maximum of the EventInfo arb_value
    AVERAGE_VALUE = 4          ///< Computes the average of the EventInfo arb_value
  };
  struct EventInfo;
  struct Event;

private:
  std::vector<RepeatingEvent> repeating_events;

public:
  size_t GetRepeatingEventTag(std::string event_name);
  void   LogEvent(size_t ev_tag,
                  EventType ev_type,
                  const std::shared_ptr<EventInfo>& ev_info);
  void   LogEvent(size_t ev_tag,
                  EventType ev_type);
  std::string PrintEventHistory(size_t ev_tag);
  double ProcessEvent(size_t ev_tag, EventOperation ev_operation);
};

//###################################################################
/** */
struct ChiLog::EventInfo
{
  std::string arb_info;
  double      arb_value = 0.0;
           EventInfo() : arb_info(std::string()) {}

  explicit EventInfo(std::string in_string) : arb_info(in_string) {}
  explicit EventInfo(double in_value) : arb_value(in_value) {}
  explicit EventInfo(std::string in_string,
                     double in_value) :
                     arb_info(in_string),
                     arb_value(in_value){}

  virtual ~EventInfo() = default;

  virtual std::string GetString()
  {
    return arb_info;
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