#include "chi_console.h"

#include "chi_runtime.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include <iostream>



//############################################################################# Run Console loop
/** Executes the loop for the console.*/
void ChiConsole::RunConsoleLoop(char*)
{

  chi_log.Log(LOG_0) << "Console loop started. "
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
      chi_log.Log(LOG_ALL) << lua_tostring(consoleState,-1);
      lua_pop(consoleState,1);
    }
  }
  chi::run_time::termination_posted = true;

  chi_log.Log(LOG_0) << "Console loop stopped successfully.";
}
