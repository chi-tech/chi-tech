#ifndef _chi_log_h
#define _chi_log_h

#include "chi_logstream.h"

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
/**Object for controlling log output. There are three levels of
 * verbosity in ChiTech: Zero(Default), One and Two. These can
 * be set on the command line via the switch -v followed by a
 * space and the number for the verbosity (0,1 or 2). The lua command
 * chiLogSetVerbosity(int_level) achieves the same.*/
class ChiLog
{
private:
  DummyStream dummy_stream;
  int verbosity;

public:
  //00
                  ChiLog();
  //01
  LogStream  Log(LOG_LVL level);
  void            SetVerbosity(int int_level);
  int             GetVerbosity();
};

#endif