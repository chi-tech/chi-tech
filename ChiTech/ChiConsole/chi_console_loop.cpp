#include "ChiConsole/chi_console.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#include <iostream>

//######################################################### Run Console loop
/** Executes the loop for the console.*/
void chi_objects::ChiConsole::RunConsoleLoop(char*) const
{
  chi::log.Log() << "Console loop started. "
                     << "Type \"exit\" to quit (or Ctl-C).";
  //======================================== Home location loop
  if (chi::mpi.location_id == 0)
  {
    while (not chi::run_time::termination_posted)
    {
      std::string console_input;

      std::cin >> console_input;

      if (console_input == std::string("exit"))
      {
        int exit_code = -1;
        MPI_Bcast(&exit_code,      //buffer
                  1, MPI_INT,      //count + type
                  0,               //root
                  MPI_COMM_WORLD); //communicator
        break;
      }

      int console_input_len = static_cast<int>(console_input.size());
      MPI_Bcast(&console_input_len, //buffer
                1, MPI_INT,         //count + type
                0,                  //root
                MPI_COMM_WORLD);    //communicator

      MPI_Bcast(console_input.data(),           //buffer
                console_input_len, MPI_CHAR,    //count + type
                0,                              //root
                MPI_COMM_WORLD);                //communicator

      if (luaL_dostring(consoleState,console_input.c_str()))
      {
        chi::log.LogAll() << lua_tostring(consoleState,-1);
        lua_pop(consoleState,1);
      }
    }
  }//if home
  //======================================== Non-home location loop
  else
  {
    while (not chi::run_time::termination_posted)
    {
      int console_input_len = -1;
      MPI_Bcast(&console_input_len,    //buffer
                1, MPI_INT,            //count + type
                0,                     //root
                MPI_COMM_WORLD);       //communicator

      if (console_input_len < 0) break;
      else
      {
        char* console_input_raw = new char[console_input_len+1];
        MPI_Bcast(console_input_raw,           //buffer
                  console_input_len, MPI_CHAR, //count + type
                  0,                           //root
                  MPI_COMM_WORLD);             //communicator
        console_input_raw[console_input_len] = '\0';

        std::string console_input(console_input_raw);

        if (luaL_dostring(consoleState,console_input.c_str()))
        {
          chi::log.LogAll() << lua_tostring(consoleState,-1);
          lua_pop(consoleState,1);
        }
      }
    }
  }

  chi::run_time::termination_posted = true;

  chi::log.Log() << "Console loop stopped successfully.";
}
