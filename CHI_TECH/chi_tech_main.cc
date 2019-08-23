#include<iostream>
#include<mpi.h>

#include"chi_tech_main.h"

void ParseArguments(int argc, char** argv);
#ifdef CHI_USEGRAPHICS
void RunInteractive(int argc, char** argv);
#endif
void RunBatch(int argc, char** argv);

/// @file


//######################################################### Program entry point
/** Program entry point.

\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.

*/
int main(int argc, char** argv)
{
  for (int k=0; k<20; k++)
    chi_global_timings[k] = 0.0;

  ParseArguments(argc,argv);
#ifdef CHI_USEGRAPHICS
  if (sim_option_withgraphics)
  {
    RunInteractive(argc,argv);
  } else
  {
    RunBatch(argc,argv);
  }
#else
  RunBatch(argc,argv);
#endif
  return 0;
}

//############################################### Small utility for parsing
/**Checks if all the characters in a
 * string is a number.*/
bool IsNumber(std::string s)
{
  for (int i = 0; i < s.length(); i++)
    if (isdigit(s[i]) == false)
      return false;

  return true;
}

//############################################### Argument parser
/**Parses input arguments.*/
void ParseArguments(int argc, char** argv)
{
  bool input_file_found = false;
  for (int i=1; i<argc; i++)
  {
    //std::cout << i;
    std::string argument(argv[i]);

    if ((argument.find('=') == std::string::npos) && (!input_file_found) )
    {
      //std::cout << "Input File: " << argument << std::endl;
      input_file_name = argument;
      input_file_found = true;
    }
    else if (argument.find('=') != std::string::npos)
    {
      luaL_dostring(chi_console.consoleState, argument.c_str());
    }
    //================================================ No-graphics option
    if ((argument.find("-b")!=std::string::npos)  )
    {
      sim_option_withgraphics = false;
    }
    //================================================ Verbosity
    if (argument.find("-v") != std::string::npos)
    {
      if ((i+1) >= argc)
      {
        std::cerr << "Invalid option used with command line argument "
                     "-v. Options are 0,1 or 2." << std::endl;
        exit(EXIT_FAILURE);
      }
      else
      {
        std::string v_option(argv[i+1]);
        try {
          int level = std::stoi(v_option);
          chi_log.SetVerbosity(level);
        }
        catch (std::invalid_argument e)
        {
          std::cerr << "Invalid option used with command line argument "
                       "-v. Options are 0,1 or 2." << std::endl;
          exit(EXIT_FAILURE);
        }
      }

    }
  }
}

#ifdef CHI_USEGRAPHICS
//========================================================= Interactive interface
/**Runs the interactive chitech engine*/
void RunInteractive(int argc, char** argv)
{
  printf("Running interactive ChiTech\n");
  chi_program_timer.Reset();
  chiwindowManager.CreateCHIWindow("ChiTech Window", false);

  if (input_file_name.size()>0)
    chi_console.ExecuteFile(input_file_name.c_str(),argc,argv);

  printf("Starting threads\n");
    omp_set_dynamic(0);
#pragma omp parallel num_threads(4)
    {
#pragma omp master
        {
            chiwindowManager.RunEventLoop();
        }
#pragma omp single nowait
        {
            chigraphics.RunGraphicsLoop(); // Graphics loop
        }
#pragma omp single nowait
        {
            chi_console.RunConsoleLoop();
        }
#pragma omp single nowait
        {
            chi_physics_handler.RunPhysicsLoop();
        }

    }
}
#endif

//############################################### Batch interface
/**Runs ChiTech in pure batch mode. Start then finish.*/
void RunBatch(int argc, char** argv)
{
  int location_id, number_processes;

  MPI_Init (&argc, &argv);                            /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &location_id);      /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &number_processes); /* get number of processes */

  chi_mpi.location_id = location_id;
  chi_mpi.process_count = number_processes;

  chi_log.Log(LOG_0)
    << "ChiTech number of arguments supplied: "
    << argc - 1;

  chi_log.Log(LOG_0)
    << "Running ChiTech in batch mode with "
    << number_processes << " processes.";

  if (argc<=1)
    chi_log.Log(LOG_0)
      << "\nUsage: exe inputfile [options values]\n"
      << "\n"
      << "     -v    Level of verbosity. Default 0. Can be either 0, 1 or 2.\n"
      << "     a=b   Executes argument as a lua string.\n\n\n";



  chi_console.PostMPIInfo(location_id, number_processes);
  chi_mpi.Initialize();

  chi_physics_handler.InitPetSc(argc,argv);

  if ((input_file_name.size()>0) )
  {
    chi_console.ExecuteFile(input_file_name.c_str(),argc,argv);
  }
  PetscFinalize();
  MPI_Finalize();
}

//############################################### Interactive interface
/**Runs the interactive chitech engine*/
void RunInteractive(int argc, char** argv)
{
  int location_id, number_processes;

  MPI_Init (&argc, &argv);                            /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &location_id);      /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &number_processes); /* get number of processes */

  chi_mpi.location_id = location_id;
  chi_mpi.process_count = number_processes;

  chi_log.Log(LOG_0)
    << "ChiTech interactive mode with " << number_processes << " processes.";

  chi_console.PostMPIInfo(location_id, number_processes);
  chi_mpi.Initialize();

  chi_physics_handler.InitPetSc(argc,argv);

  if ((input_file_name.size()>0) )
  {
    chi_console.ExecuteFile(input_file_name.c_str(),argc,argv);
  }



  PetscFinalize();
  MPI_Finalize();
}