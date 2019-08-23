#include"chi_timer.h"


//############################################################################# Default constr
/** Default constructor.*/
CHI_TIMER::CHI_TIMER()
{
    #ifdef UNIX_ENV
      clock_gettime(CLOCK_MONOTONIC,&this->startTime);
    #else
      this->CounterStart = (__int64)0.0;
	    this->counterTime = 0.0;
	    this->PCFreq = 1000.0;

      QueryPerformanceFrequency(&this->LargeInt);
	    this->PCFreq = double(this->LargeInt.QuadPart)/1000.0;
    #endif
}
