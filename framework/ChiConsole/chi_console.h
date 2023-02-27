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
namespace chi_objects
{
/** Class for handling the console and scripting.*/
class ChiConsole
{
private:
  lua_State*							consoleState;    ///< Pointer to lua console state

  std::vector<std::string> command_buffer_; ///< Buffer of commands to execute
  static ChiConsole       instance_;
  //00
              ChiConsole() noexcept;
public:
  static ChiConsole& GetInstance() noexcept
  {return instance_;}

  lua_State*& GetConsoleState() {return consoleState;}
  std::vector<std::string>& GetCommandBuffer() {return command_buffer_;}

  //01 Loop
  void        RunConsoleLoop(char* fileName=nullptr) const;
  //02 Utilities
  int         ExecuteFile(const std::string& fileName,int argc, char** argv) const;
  void        PostMPIInfo(int location_id, int number_of_processes) const;
  void        RegisterFunction(const std::string& string_name,
                               lua_CFunction function);
  //03
  void        FlushConsole();
  //05 Memory
  static CSTMemory  GetMemoryUsage();
  static double      GetMemoryUsageInMB();
};
}//namespace chi_objects


#endif
