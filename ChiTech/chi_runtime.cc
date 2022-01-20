/** @file Runtime file*/
#include "chi_runtime.h"

#include "ChiConsole/chi_console.h"
#include "ChiMath/chi_math.h"
#include "ChiPhysics/chi_physics.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "chi_mpi.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include <iostream>

//=============================================== Global variables
chi_math::UnknownManager ChiMath::UNITARY_UNKNOWN_MANAGER;

ChiConsole  ChiConsole::instance;
ChiMath     ChiMath::instance;
ChiMPI      ChiMPI::instance;
ChiLog      ChiLog::instance;
ChiPhysics  ChiPhysics::instance;



ChiConsole&  chi_console = ChiConsole::GetInstance();
ChiMath&     chi_math_handler = ChiMath::GetInstance();
ChiMPI&      chi_mpi = ChiMPI::GetInstance();
ChiLog&      chi_log = ChiLog::GetInstance();
ChiPhysics&  chi_physics_handler = ChiPhysics::GetInstance();

ChiTimer    chi_program_timer;

/** Global stack of handlers */
std::vector<chi_mesh::MeshHandler*>  ChiTech::meshhandler_stack;
int                                  ChiTech::current_mesh_handler=-1;
bool                                 ChiTech::termination_posted = false;
std::string                          ChiTech::input_file_name;
bool                                 ChiTech::sim_option_interactive = true;
bool                                 ChiTech::allow_petsc_error_handler = false;




//############################################### Argument parser
/**Parses input arguments.
\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.
 */
void ChiTech::ParseArguments(int argc, char** argv)
{
  bool input_file_found = false;
  for (int i=1; i<argc; i++)
  {
    std::string argument(argv[i]);

    chi_log.Log(LOG_0) << "Parsing argument " << i << " " << argument;

    if (argument.find("-h")!=std::string::npos)
    {
      chi_log.Log(LOG_0)
        << "\nUsage: exe inputfile [options values]\n"
        << "\n"
        << "     -v                         Level of verbosity. Default 0. Can be either 0, 1 or 2.\n"
        << "     a=b                        Executes argument as a lua string. i.e. x=2 or y=[[\"string\"]]\n"
        << "     -allow_petsc_error_handler Allow petsc error handler.\n\n\n";

      chi_log.Log(LOG_0) << "PETSc options:";
      ChiTech::termination_posted = true;
    }
    else if (argument.find("-allow_petsc_error_handler")!=std::string::npos)
    {
      ChiTech::allow_petsc_error_handler = true;
    }
    //================================================ No-graphics option
    else if (argument.find("-b")!=std::string::npos)
    {
      ChiTech::sim_option_interactive = false;
    }//-b
    //================================================ Verbosity
    else if (argument.find("-v") != std::string::npos)
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
        catch (const std::invalid_argument& e)
        {
          std::cerr << "Invalid option used with command line argument "
                       "-v. Options are 0,1 or 2." << std::endl;
          exit(EXIT_FAILURE);
        }
      }

    }//-v
    else if ((argument.find('=') == std::string::npos) && (!input_file_found) )
    {
      ChiTech::input_file_name = argument;
      input_file_found = true;
      ChiTech::sim_option_interactive = false;
    }//no =
    else if (argument.find('=') != std::string::npos)
    {
      chi_console.command_buffer.push_back(argument);
    }//=

  }//for argument
}

//############################################### Initialize ChiTech
/**Initializes all necessary items for ChiTech.
\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.
 */
int ChiTech::Initialize(int argc, char** argv)
{
  int location_id = 0, number_processes = 1;

  MPI_Init (&argc, &argv);                           /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &location_id);      /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &number_processes); /* get number of processes */

  chi_mpi.SetLocationID(location_id);
  chi_mpi.SetProcessCount(number_processes);

  chi_console.PostMPIInfo(location_id, number_processes);

  ParseArguments(argc, argv);

  chi_physics_handler.InitPetSc(argc,argv);

  return 0;
}

//############################################### Finalize ChiTech
/**Finalizes ChiTech.
 * */
void ChiTech::Finalize()
{
  PetscFinalize();
  MPI_Finalize();
}

//############################################### Interactive interface
/**Runs the interactive chitech engine*/
int ChiTech::RunInteractive(int argc, char** argv)
{
  chi_log.Log(LOG_0)
    << ChiTimer::GetLocalDateTimeString()
    << " Running ChiTech in interactive-mode with "
    << chi_mpi.process_count << " processes.";

  chi_log.Log(LOG_0)
    << "ChiTech number of arguments supplied: "
    << argc - 1;

  chi_console.FlushConsole();

  if ( not ChiTech::input_file_name.empty() )
    chi_console.ExecuteFile(ChiTech::input_file_name.c_str(), argc, argv);

  chi_console.RunConsoleLoop();

  chi_log.Log(LOG_0)
    << "Final program time " << chi_program_timer.GetTimeString();

  chi_log.Log(LOG_0)
    << ChiTimer::GetLocalDateTimeString()
    << " ChiTech finished execution.";

  return 0;
}

//############################################### Batch interface
/**Runs ChiTech in pure batch mode. Start then finish.*/
int ChiTech::RunBatch(int argc, char** argv)
{
  chi_log.Log(LOG_0)
    << ChiTimer::GetLocalDateTimeString()
    << " Running ChiTech in batch-mode with "
    << chi_mpi.process_count << " processes.";

  chi_log.Log(LOG_0)
    << "ChiTech number of arguments supplied: "
    << argc - 1;

  if (argc<=1)
    chi_log.Log(LOG_0)
      << "\nUsage: exe inputfile [options values]\n"
      << "\n"
      << "     -v                         Level of verbosity. Default 0. Can be either 0, 1 or 2.\n"
      << "     a=b                        Executes argument as a lua string. i.e. x=2 or y=[[\"string\"]]\n"
      << "     -allow_petsc_error_handler Allow petsc error handler.\n\n\n";

  chi_console.FlushConsole();

  int error_code = 0;
  if ( not ChiTech::input_file_name.empty() )
    error_code = chi_console.ExecuteFile(ChiTech::input_file_name.c_str(), argc, argv);

  chi_log.Log(LOG_0)
    << "Final program time " << chi_program_timer.GetTimeString();

  chi_log.Log(LOG_0)
    << ChiTimer::GetLocalDateTimeString()
    << " ChiTech finished execution of " << ChiTech::input_file_name;

  return error_code;
}


