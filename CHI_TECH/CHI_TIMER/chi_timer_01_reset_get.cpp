#include"chi_timer.h"

#include <math.h>

//################################################################### Reset
/** Resets the timer to zero.*/
void CHI_TIMER::Reset()
{
    #ifdef UNIX_ENV
	     clock_gettime(CLOCK_MONOTONIC,&this->startTime);
    #else
	     QueryPerformanceCounter(&this->LargeInt);
	     CounterStart = this->LargeInt.QuadPart;
    #endif
}

//################################################################### Get time
/** Gets the current timer value in milliseconds.*/
double CHI_TIMER::GetTime()
{
    #ifdef UNIX_ENV
	     timespec  newTime;
	     clock_gettime(CLOCK_MONOTONIC,&newTime);
	     timespec diff=this->GetDifference(this->startTime,newTime);
	     double diffTime=diff.tv_sec*1000.0+diff.tv_nsec/1000000.0;
    #else
	     QueryPerformanceFrequency(&this->LargeInt);
		   this->PCFreq = double(this->LargeInt.QuadPart)/1000.0;
       QueryPerformanceCounter(&this->LargeInt);
	     this->counterTime = double(this->LargeInt.QuadPart - this->CounterStart)/this->PCFreq;
       double diffTime=this->counterTime;
    #endif

	return diffTime;
}

//################################################################### Get string
/**Obtains a stringstream in the format of hh:mm::ss.
 *
 * \param time_milli Time in milliseconds.
 *
 * */
std::string CHI_TIMER::GetTimeString()
{
  double time_sec = this->GetTime()/1000.0;
  int    hours    = std::floor(time_sec/60/60);
  int    minutes  = std::floor((time_sec-60*60*hours)/60);
  int    seconds  = time_sec - 3600*hours - 60*minutes;

  char buff[100];
  sprintf(buff,"%02d:%02d:%02d",hours,minutes,seconds);

  return std::string(buff);
}


#ifdef UNIX_ENV
//################################################################### Get Difference
/** Determines the difference between two time calls.*/
timespec CHI_TIMER::GetDifference(timespec start, timespec end)
{

	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;

    //return 0.0;

}

#endif
