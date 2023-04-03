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
#include <map>

class chi;

/**Small utility macro for joining two words.*/
#define ChiConsoleJoinWordsA(x,y) x##y
/**IDK why this is needed. Seems like counter doesnt work properly without it*/
#define ChiConsoleJoinWordsB(x,y) ChiConsoleJoinWordsA(x,y)

/**Macro for registering a lua_CFunction within the ChiConsole
 * singleton, with the function being in the global namespace. Example:
 * \code
 * ChiConsoleRegisterLuaFunction("chiSolverInitialize");
 * \endcode
 *
 * \note Remember to include the header "ChiConsole/chi_console.h"*/
#define ChiConsoleRegisterLuaFunction(func_name) \
  static char ChiConsoleJoinWordsB(               \
  unique_var_name_luacfunc_##func_name##_, __COUNTER__) = \
    chi_objects::ChiConsole::AddFunctionToRegistry(#func_name, func_name)

//############################################################################# CLASS DEF
namespace chi_objects
{

/** Class for handling the console and scripting.*/
class ChiConsole
{
private:
  struct LuaFunctionRegistryEntry
  {
    lua_CFunction function_ptr;
    std::string   function_raw_name;
  };
private:
  lua_State*							console_state_;    ///< Pointer to lua console state

  std::vector<std::string> command_buffer_; ///< Buffer of commands to execute
  static ChiConsole       instance_;

  std::map<std::string, LuaFunctionRegistryEntry> lua_function_registry_;
  //00
  ChiConsole() noexcept;
private:
  friend class ::chi;
  void LoadRegisteredLuaItems();
public:
  static ChiConsole& GetInstance() noexcept;

  lua_State*& GetConsoleState() {return console_state_;}
  std::vector<std::string>& GetCommandBuffer() {return command_buffer_;}

  //01 Loop
  void        RunConsoleLoop(char* fileName=nullptr) const;
  //02 Utilities
  int         ExecuteFile(const std::string& fileName,int argc, char** argv) const;
  void        PostMPIInfo(int location_id, int number_of_processes) const;
  static char AddFunctionToRegistry(const std::string& raw_name_in_lua,
                                    lua_CFunction function_ptr);
  //03
  void        FlushConsole();
  //05 Memory
  static CSTMemory  GetMemoryUsage();
  static double      GetMemoryUsageInMB();
};
}//namespace chi_objects


#endif
