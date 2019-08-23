#ifndef CHI_CONSOLE_H
#define CHI_CONSOLE_H

#include "../CHI_VECTOR/chi_vector.h"
extern "C"
{
#include<lua.h>
#include<lualib.h>
#include<lauxlib.h>
}
#include "chi_console_structs.h"


//template class CHI_VECTOR<CST_EVENT>;

//############################################################################# CLASS DEF
/** Class for handling the console and scripting.*/
class CHI_CONSOLE
{
    public:
	lua_State*							consoleState;             	///< Pointer to lua console state
	CHI_VECTOR<CST_EVENT>		eventsStack;         		///< A vector of all events occurring
	bool                    exitLoop;
	bool                    runDeveloper;
	long                    currentSize;
	long                    previousSize;
	char                    buffer[2000];
	int                     numberOfLines;
	int                     xSize;
    public:
	//00
						  CHI_CONSOLE();
  //01 Loop
  void        RunConsoleLoop(char* fileName=NULL);
  //02 Utilities
  void        InitializeLuaEvent(const char* eventTitle);
  void        ExecuteFile(const char* fileName,int argc, char** argv);
  void				PostEventToConsole(CST_EVENT* inputEvent);
  void        PostMPIInfo(int location_id, int number_of_processes);
  //03
  void        flushConsole();
  //04 Console Info
  int         GetNumCharsInConsoleBuffer();
  void        CopyConsole(char* destination,int lineNumber=0,int xSize=80);
  //05 Memory
  CST_MEMORY  GetMemoryUsage();
  double      GetMemoryUsageInMB();
  double      GetMemoryUsageInBytes();

};



#endif
