#ifndef CHI_CONSOLE_H
#define CHI_CONSOLE_H

extern "C"
{
#include<lua.h>
#include<lualib.h>
#include<lauxlib.h>
}
#include "chi_console_structs.h"


//template class CHI_VECTOR<CSTEvent>;

//############################################################################# CLASS DEF
/** Class for handling the console and scripting.*/
class ChiConsole
{
    public:
	lua_State*							consoleState;             	///< Pointer to lua console state
	bool                    exitLoop;
	bool                    runDeveloper;
	long                    currentSize;
	long                    previousSize;
	char                    buffer[2000];
	int                     numberOfLines;
	int                     xSize;
    public:
	//00
						  ChiConsole();
  //01 Loop
  void        RunConsoleLoop(char* fileName=NULL);
  //02 Utilities
  void        InitializeLuaEvent(const char* eventTitle);
  void        ExecuteFile(const char* fileName,int argc, char** argv);
  void				PostEventToConsole(CSTEvent* inputEvent);
  void        PostMPIInfo(int location_id, int number_of_processes);
  //03
  void        flushConsole();
  //04 Console Info
  int         GetNumCharsInConsoleBuffer();
  void        CopyConsole(char* destination,int lineNumber=0,int xSize=80);
  //05 Memory
  CSTMemory  GetMemoryUsage();
  double      GetMemoryUsageInMB();
  double      GetMemoryUsageInBytes();

};



#endif
