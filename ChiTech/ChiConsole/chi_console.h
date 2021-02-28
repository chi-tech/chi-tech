#ifndef CHI_CONSOLE_H
#define CHI_CONSOLE_H

extern "C"
{
#include<lua.h>
#include<lualib.h>
#include<lauxlib.h>
}
#include "chi_console_structs.h"

#include <vector>
#include <string>

//############################################################################# CLASS DEF
/** Class for handling the console and scripting.*/
class ChiConsole
{
public:
	lua_State*							consoleState;             	///< Pointer to lua console state
	bool                    exit_loop=false;
	bool                    runDeveloper=false;
	long                    currentSize=0;
	long                    previousSize=0;
	char                    buffer[2000]={'\0'};
	int                     numberOfLines=0;
	int                     xSize=80;

	std::vector<std::string> command_buffer;

private:
  static ChiConsole       instance;
	//00
						  ChiConsole() noexcept;
public:
	static ChiConsole& GetInstance() noexcept
  {return instance;}

  //01 Loop
  void        RunConsoleLoop(char* fileName=nullptr);
  //02 Utilities
  void        InitializeLuaEvent(const char* eventTitle);
  int         ExecuteFile(const char* fileName,int argc, char** argv);
  void				PostEventToConsole(CSTEvent* inputEvent);
  void        PostMPIInfo(int location_id, int number_of_processes);
  //03
  void        FlushConsole();
  //04 Console Info
  int         GetNumCharsInConsoleBuffer();
  void        CopyConsole(char* destination,int lineNumber=0,int xSize=80);
  //05 Memory
  CSTMemory  GetMemoryUsage();
  double      GetMemoryUsageInMB();
  double      GetMemoryUsageInBytes();

};



#endif
