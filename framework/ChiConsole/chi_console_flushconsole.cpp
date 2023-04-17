#include "ChiConsole/chi_console.h"
#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/** This function sends the commands contained in the command_buffer to
the lua state from where it is executed. These could be commands passed
 via the command line or loaded elsewhere.*/
void chi_objects::ChiConsole::FlushConsole()
{
  try
  {
    for (auto& command : command_buffer_)
    {
      bool error = luaL_dostring(console_state_, command.c_str());
      if (error)
      {
        chi::log.LogAll() << lua_tostring(console_state_, -1);
        lua_pop(console_state_, 1);
      }
    }
  }
  catch(const std::exception& e)
  {
    chi::log.LogAllError() << e.what();
    chi::Exit(EXIT_FAILURE);
  }
}
