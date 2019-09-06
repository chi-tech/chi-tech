#ifndef CHI_TIMER_H
#define CHI_TIMER_H

#include <string>

#ifdef UNIX_ENV
	#include<time.h>
#else
	#include<time.h>
	#include<Windows.h>
	#include<stddef.h>
	#include<iostream>
	#include<fstream>
#endif

#ifdef UNIX_ENV
//############################################################################# CLASS DEF
/** Timer object.*/
class CHI_TIMER
{
public:
	timespec  	startTime;
public:
	//00
				CHI_TIMER();
	//01
	void   		Reset();
	double 		GetTime();
	std::string GetTimeString();
	timespec 	GetDifference(timespec start, timespec end);
};
#else

//######################################################### CLASS DEFINITION
/** Timer object.*/
class CHI_TIMER
{
public:
	double			PCFreq;								///< Computer tick frequency in Hz
	double			counterTime;
	__int64			CounterStart;						///< Start ticks registered by timer
	LARGE_INTEGER	LargeInt;							///< Type declaration of large integer

public:
    //00
					CHI_TIMER();
    //01
	void			Reset();
	std::string GetTimeString();
	double			GetTime();
};
#endif

#endif
