#ifndef CHI_TIMER_H
#define CHI_TIMER_H

#include <string>
#include <chrono>

//############################################################################# CLASS DEF
/** Timer object.*/
class ChiTimer
{
public:
  std::chrono::steady_clock::time_point startTime;

public:
	//00
				      ChiTimer() noexcept;
	//01
	void   		  Reset();
	double 		  GetTime();
	std::string GetTimeString();
	std::string GetLocalDateTimeString();
};

#endif

