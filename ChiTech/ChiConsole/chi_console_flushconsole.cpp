#include "chi_console.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&   		chi_physics_handler;

#include "ChiTimer/chi_timer.h"
extern ChiTimer        chi_program_timer;

//###################################################################
/** This function sends the commands contained in the command_buffer to
the lua state from where it is executed. These could be commands passed
 via the command line or loaded elsewhere.*/
void ChiConsole::FlushConsole()
{
  for (auto& command : command_buffer)
    luaL_dostring(consoleState, command.c_str());
}
