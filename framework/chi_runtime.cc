/** @file Runtime file*/
#include "chi_runtime.h"
#include "chi_configuration.h"

#include "ChiConsole/chi_console.h"
#include "ChiMath/chi_math.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "ChiPhysics/chi_physics_namespace.h"

#include "ChiObject/object_maker.h"

#include "chi_mpi.h"
#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

#include <iostream>

#ifndef NDEBUG
#include <unistd.h>
#endif

//=============================================== Global variables
chi::ChiConsole& Chi::console = chi::ChiConsole::GetInstance();
chi::ChiLog& Chi::log = chi::ChiLog::GetInstance();
chi::MPI_Info& Chi::mpi = chi::MPI_Info::GetInstance();
chi::ChiTimer Chi::program_timer;

/** Global stack of handlers */
std::vector<chi_mesh::MeshHandlerPtr> Chi::meshhandler_stack;
int Chi::current_mesh_handler = -1;

std::vector<chi_mesh::SurfaceMeshPtr> Chi::surface_mesh_stack;
std::vector<chi_mesh::FFInterpPtr> Chi::field_func_interpolation_stack;
std::vector<chi_mesh::UnpartMeshPtr> Chi::unpartitionedmesh_stack;

std::vector<chi_physics::MaterialPtr> Chi::material_stack;
std::vector<chi_physics::MultiGroupXSPtr> Chi::multigroup_xs_stack;
std::vector<chi_physics::FieldFunctionPtr> Chi::field_function_stack;

std::vector<chi_math::QuadraturePtr> Chi::quadrature_stack;
std::vector<chi_math::AngularQuadraturePtr> Chi::angular_quadrature_stack;

std::vector<ChiObjectPtr> Chi::object_stack;
std::vector<chi_math::SpatialDiscretizationPtr> Chi::sdm_stack;

//================================ run_time quantities
bool Chi::run_time::termination_posted_ = false;
std::string Chi::run_time::input_file_name_;
bool Chi::run_time::sim_option_interactive_ = true;
bool Chi::run_time::allow_petsc_error_handler_ = false;
bool Chi::run_time::supress_beg_end_timelog_ = false;
bool Chi::run_time::suppress_color_ = false;
bool Chi::run_time::dump_registry_ = false;

const std::string Chi::run_time::command_line_help_string_ =
  "\nUsage: exe inputfile [options values]\n"
  "\n"
  "     -v                          Level of verbosity. Default 0.\n"
  "                                 Can be either 0, 1 or 2.\n"
  "     a=b                         Executes argument as a lua string. "
  "i.e. x=2 or y=[[\"string\"]]\n"
  "     --allow_petsc_error_handler Allow petsc error handler.\n"
  "     --supress_beg_end_timelog   Suppress time logs at the \n"
  "                                 beginning and end of execution.\n"
  "     --suppress_color            Suppresses the printing of color.\n"
  "                                 useful for unit tests requiring a diff.\n"
  "     --dump-object-registry      Dumps the object registry.\n"
  "\n\n\n";

// ############################################### Argument parser
/**Parses input arguments.
\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.
 */
void Chi::run_time::ParseArguments(int argc, char** argv)
{
  bool input_file_found = false;
  for (int i = 1; i < argc; i++)
  {
    std::string argument(argv[i]);

    Chi::log.Log() << "Parsing argument " << i << " " << argument;

    if (argument.find("-h") != std::string::npos or
        argument.find("--help") != std::string::npos)
    {
      Chi::log.Log() << Chi::run_time::command_line_help_string_;
      Chi::run_time::termination_posted_ = true;
    }
    else if (argument.find("--supress_beg_end_timelog") != std::string::npos)
    {
      Chi::run_time::supress_beg_end_timelog_ = true;
    }
    else if (argument.find("--allow_petsc_error_handler") != std::string::npos)
    {
      Chi::run_time::allow_petsc_error_handler_ = true;
    }
    else if (argument.find("--suppress_color") != std::string::npos)
    {
      Chi::run_time::suppress_color_ = true;
    }
    else if (argument.find("--dump-object-registry") != std::string::npos)
    {
      Chi::run_time::dump_registry_ = true;
      Chi::run_time::termination_posted_ = true;
    }
    //================================================ No-graphics option
    else if (argument.find("-b") != std::string::npos)
    {
      Chi::run_time::sim_option_interactive_ = false;
    } //-b
    //================================================ Verbosity
    else if (argument.find("-v") != std::string::npos)
    {
      if ((i + 1) >= argc)
      {
        std::cerr << "Invalid option used with command line argument "
                     "-v. Options are 0,1 or 2."
                  << std::endl;
        Chi::Exit(EXIT_FAILURE);
      }
      else
      {
        std::string v_option(argv[i + 1]);
        try
        {
          int level = std::stoi(v_option);
          Chi::log.SetVerbosity(level);
        }
        catch (const std::invalid_argument& e)
        {
          std::cerr << "Invalid option used with command line argument "
                       "-v. Options are 0,1 or 2."
                    << std::endl;
          Chi::Exit(EXIT_FAILURE);
        }
      }

    } //-v
    else if ((argument.find('=') == std::string::npos) and (!input_file_found))
    {
      Chi::run_time::input_file_name_ = argument;
      input_file_found = true;
      Chi::run_time::sim_option_interactive_ = false;
    } // no =
    else if (argument.find('=') != std::string::npos)
    {
      Chi::console.GetCommandBuffer().push_back(argument);
    } //=

  } // for argument

  if (Chi::run_time::dump_registry_)
  {
    ChiObjectMaker::GetInstance().DumpRegister();
    Chi::console.DumpRegister();
  }

}

// ############################################### Initialize ChiTech
/**Initializes all necessary items for ChiTech.
\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.
 */
int Chi::Initialize(int argc, char** argv)
{
  int location_id = 0, number_processes = 1;

  MPI_Init(&argc, &argv);                           /* starts MPI */
  MPI_Comm_rank(MPI_COMM_WORLD, &location_id);      /* get cur process id */
  MPI_Comm_size(MPI_COMM_WORLD, &number_processes); /* get num of processes */

  mpi.SetLocationID(location_id);
  mpi.SetProcessCount(number_processes);

  Chi::console.LoadRegisteredLuaItems();
  Chi::console.PostMPIInfo(location_id,       number_processes);

  run_time::ParseArguments(argc, argv);

  run_time::InitPetSc(argc, argv);

  return 0;
}

/**Initializes PetSc for use by all entities.*/
int Chi::run_time::InitPetSc(int argc, char** argv)
{
  PetscErrorCode ierr;
  PetscMPIInt size;

  PetscOptionsInsertString(nullptr, "-error_output_stderr");
  if (not Chi::run_time::allow_petsc_error_handler_)
    PetscOptionsInsertString(nullptr, "-no_signal_handler");

  ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
  if (ierr) return ierr;

  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
  CHKERRQ(ierr);

  return 0;
}

// ############################################### Finalize ChiTech
/**Finalizes ChiTech.
 * */
void Chi::Finalize()
{
  meshhandler_stack.clear();

  surface_mesh_stack.clear();
  object_stack.clear();
  field_func_interpolation_stack.clear();
  unpartitionedmesh_stack.clear();

  object_stack.clear();
  material_stack.clear();
  multigroup_xs_stack.clear();

  PetscFinalize();
  MPI_Finalize();
}

// ############################################### Interactive interface
/**Runs the interactive chitech engine*/
int Chi::RunInteractive(int argc, char** argv)
{
  if (not Chi::run_time::supress_beg_end_timelog_)
  {
    Chi::log.Log() << chi::ChiTimer::GetLocalDateTimeString()
                   << " Running ChiTech in interactive-mode with "
                   << Chi::mpi.process_count << " processes.";

    Chi::log.Log() << "ChiTech version " << GetVersionStr();
  }

  Chi::log.Log() << "ChiTech number of arguments supplied: " << argc - 1;

  Chi::log.LogAll();

  Chi::console.FlushConsole();

  const auto& input_fname = Chi::run_time::input_file_name_;

  if (not input_fname.empty())
  {
    try
    {
      Chi::console.ExecuteFile(input_fname, argc, argv);
    }
    catch (const std::exception& excp)
    {
      Chi::log.LogAllError() << excp.what();
      // No quitting if file execution fails
    }
  }

  Chi::console.RunConsoleLoop();

  if (not Chi::run_time::supress_beg_end_timelog_)
  {
    Chi::log.Log() << "Final program time " << program_timer.GetTimeString();
    Chi::log.Log() << chi::ChiTimer::GetLocalDateTimeString()
                   << " ChiTech finished execution.";
  }

  return 0;
}

// ############################################### Batch interface
/**Runs ChiTech in pure batch mode. Start then finish.*/
int Chi::RunBatch(int argc, char** argv)
{
  if (not Chi::run_time::supress_beg_end_timelog_)
  {
    Chi::log.Log() << chi::ChiTimer::GetLocalDateTimeString()
                   << " Running ChiTech in batch-mode with "
                   << Chi::mpi.process_count << " processes.";

    Chi::log.Log() << "ChiTech version " << GetVersionStr();
  }

  Chi::log.Log() << "ChiTech number of arguments supplied: " << argc - 1;

  if (argc <= 1) Chi::log.Log() << Chi::run_time::command_line_help_string_;
  Chi::console.FlushConsole();

#ifndef NDEBUG
  Chi::log.Log() << "Waiting...";
  if (Chi::mpi.location_id == 0)
    for (int k = 0; k < 2; ++k)
    {
      usleep(1000000);
      Chi::log.Log() << k;
    }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  const auto& input_fname = Chi::run_time::input_file_name_;
  int error_code = 0;

  if ((not input_fname.empty()) and (not Chi::run_time::termination_posted_))
  {
    try
    {
      error_code = Chi::console.ExecuteFile(input_fname, argc, argv);
    }
    catch (const std::exception& excp)
    {
      Chi::log.LogAllError() << excp.what();
      Chi::Exit(EXIT_FAILURE);
    }
  }

  if (not Chi::run_time::supress_beg_end_timelog_)
  {
    Chi::log.Log() << "\nFinal program time " << program_timer.GetTimeString();
    Chi::log.Log() << chi::ChiTimer::GetLocalDateTimeString()
                   << " ChiTech finished execution of "
                   << Chi::run_time::input_file_name_;
  }

  return error_code;
}

// ###################################################################
/** Exits the program appropriately.*/
void Chi::Exit(int error_code) { MPI_Abort(MPI_COMM_WORLD, error_code); }

// ###################################################################
/** Gets the ChiTech-version string.*/
std::string Chi::GetVersionStr() { return PROJECT_VERSION; }
