#include"chi_timer.h"

#include <cmath>
#include <ctime>

//################################################################### Default constr
/** Default constructor.*/
chi_objects::ChiTimer::ChiTimer() noexcept
{
  startTime = std::chrono::steady_clock::now();
}

//################################################################### Reset
/** Resets the timer to zero.*/
void chi_objects::ChiTimer::Reset()
{
  startTime = std::chrono::steady_clock::now();
}

//################################################################### Get time
/** Gets the current timer value in milliseconds.*/
double chi_objects::ChiTimer::GetTime() const
{
  using namespace std::chrono;

  steady_clock::time_point newTime = std::chrono::steady_clock::now();
  duration<double> time_span =
    duration_cast<duration<double>>(newTime - startTime);

  return time_span.count()*1000.0;
}

//################################################################### Get string
/**Obtains a string in the format of hh:mm::ss.
 *
 * */
std::string chi_objects::ChiTimer::GetTimeString() const
{
  double time_sec = this->GetTime()/1000.0;
  int    hours    = std::floor(time_sec/60/60);
  int    minutes  = std::floor((time_sec-60*60*hours)/60);
  int    seconds  = (int)time_sec - 3600*hours - 60*minutes;

  char buff[100];
  snprintf(buff,100,"%02d:%02d:%02d",hours,minutes,seconds);

  return {buff};
}

//################################################################### Get date
/**Obtains a string in the format YYYY-MM-DD hh:mm:ss
 *
 * */
std::string chi_objects::ChiTimer::GetLocalDateTimeString()
{
  using namespace std::chrono;
  std::time_t now = system_clock::to_time_t(system_clock::now());

  char s[30];
  size_t end = std::strftime(s, 30, "%Y-%m-%d %H:%M:%S", std::localtime(&now));
  s[29] = '\0';
  if (end < 30)  s[end]='\0';
  return s;
}
