#include "ChiConsole/chi_console.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <iostream>



//############################################################################# Run Console loop
/** Executes the loop for the console.*/
void chi_objects::ChiConsole::RunConsoleLoop(char*)
{

  chi::log.Log() << "Console loop started. "
                     << "Type \"exit\" to quit (or Ctl-C).";
  exit_loop = false;

  while ((!exit_loop) and (!chi::run_time::termination_posted))
  {
    std::string console_input;

    std::cin >> console_input;

    if (console_input == std::string("exit"))
    {
      exit_loop=true;
      break;
    }

    if (luaL_dostring(consoleState,console_input.c_str()))
    {
      chi::log.LogAll() << lua_tostring(consoleState,-1);
      lua_pop(consoleState,1);
    }
  }
  chi::run_time::termination_posted = true;

  chi::log.Log() << "Console loop stopped successfully.";
}
