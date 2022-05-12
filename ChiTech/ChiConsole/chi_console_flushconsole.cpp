#include "ChiConsole/chi_console.h"

//###################################################################
/** This function sends the commands contained in the command_buffer to
the lua state from where it is executed. These could be commands passed
 via the command line or loaded elsewhere.*/
void chi_objects::ChiConsole::FlushConsole()
{
  for (auto& command : command_buffer)
    luaL_dostring(consoleState, command.c_str());
}
