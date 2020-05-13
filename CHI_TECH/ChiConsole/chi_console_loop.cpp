#include "chi_console.h"
#include <iostream>

extern bool chi_termination_posted;

#include <chi_log.h>
extern ChiLog& chi_log;


//############################################################################# Run Console loop
/** Executes the loop for the console.*/
void ChiConsole::RunConsoleLoop(char* fileName)
{

  chi_log.Log(LOG_0) << "Console loop started.";
  exit_loop = false;

  while ((!exit_loop) and (!chi_termination_posted))
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
  chi_termination_posted = true;

  chi_log.Log(LOG_0) << "Console loop stopped successfully.";
}
