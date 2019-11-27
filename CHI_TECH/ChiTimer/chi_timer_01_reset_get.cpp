#include"chi_timer.h"

#include <cmath>

//################################################################### Reset
/** Resets the timer to zero.*/
void ChiTimer::Reset()
{
  startTime = std::chrono::steady_clock::now();
}

//################################################################### Get time
/** Gets the current timer value in milliseconds.*/
double ChiTimer::GetTime()
{
  using namespace std::chrono;

  steady_clock::time_point newTime = std::chrono::steady_clock::now();
  duration<double> time_span =
    duration_cast<duration<double>>(newTime - startTime);

	return time_span.count()*1000.0;
}

//################################################################### Get string
/**Obtains a stringstream in the format of hh:mm::ss.
 *
 * \param time_milli Time in milliseconds.
 *
 * */
std::string ChiTimer::GetTimeString()
{
  double time_sec = this->GetTime()/1000.0;
  int    hours    = std::floor(time_sec/60/60);
  int    minutes  = std::floor((time_sec-60*60*hours)/60);
  int    seconds  = (int)time_sec - 3600*hours - 60*minutes;

  char buff[100];
  sprintf(buff,"%02d:%02d:%02d",hours,minutes,seconds);

  return std::string(buff);
}


