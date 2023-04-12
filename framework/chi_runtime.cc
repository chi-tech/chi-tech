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
chi_objects::ChiConsole& chi::console = chi_objects::ChiConsole::GetInstance();
chi_objects::ChiLog& chi::log = chi_objects::ChiLog::GetInstance();
chi_objects::MPI_Info& chi::mpi = chi_objects::MPI_Info::GetInstance();
chi_objects::ChiTimer chi::program_timer;

/** Global stack of handlers */
std::vector<chi_mesh::MeshHandlerPtr> chi::meshhandler_stack;
int chi::current_mesh_handler = -1;

std::vector<chi_mesh::SurfaceMeshPtr> chi::surface_mesh_stack;
std::vector<chi_mesh::LogicalVolumePtr> chi::logicvolume_stack;
std::vector<chi_mesh::FFInterpPtr> chi::field_func_interpolation_stack;
std::vector<chi_mesh::UnpartMeshPtr> chi::unpartitionedmesh_stack;

std::vector<chi_physics::MaterialPtr> chi::material_stack;
std::vector<chi_physics::MultiGroupXSPtr> chi::multigroup_xs_stack;
std::vector<chi_physics::FieldFunctionPtr> chi::field_function_stack;

std::vector<chi_math::QuadraturePtr> chi::quadrature_stack;
std::vector<chi_math::AngularQuadraturePtr> chi::angular_quadrature_stack;

std::vector<ChiObjectPtr> chi::object_stack;

//================================ run_time quantities
bool chi::run_time::termination_posted_ = false;
std::string chi::run_time::input_file_name_;
bool chi::run_time::sim_option_interactive_ = true;
bool chi::run_time::allow_petsc_error_handler_ = false;
bool chi::run_time::supress_beg_end_timelog_ = false;
bool chi::run_time::suppress_color_ = false;
bool chi::run_time::dump_registry_ = false;

const std::string chi::run_time::command_line_help_string_ =
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
void chi::run_time::ParseArguments(int argc, char** argv)
{
  bool input_file_found = false;
  for (int i = 1; i < argc; i++)
  {
    std::string argument(argv[i]);

    chi::log.Log() << "Parsing argument " << i << " " << argument;

    if (argument.find("-h") != std::string::npos or
        argument.find("--help") != std::string::npos)
    {
      chi::log.Log() << chi::run_time::command_line_help_string_;
      chi::run_time::termination_posted_ = true;
    }
    else if (argument.find("--supress_beg_end_timelog") != std::string::npos)
    {
      chi::run_time::supress_beg_end_timelog_ = true;
    }
    else if (argument.find("--allow_petsc_error_handler") != std::string::npos)
    {
      chi::run_time::allow_petsc_error_handler_ = true;
    }
    else if (argument.find("--suppress_color") != std::string::npos)
    {
      chi::run_time::suppress_color_ = true;
    }
    else if (argument.find("--dump-object-registry") != std::string::npos)
    {
      chi::run_time::dump_registry_ = true;
      chi::run_time::termination_posted_ = true;
    }
    //================================================ No-graphics option
    else if (argument.find("-b") != std::string::npos)
    {
      chi::run_time::sim_option_interactive_ = false;
    } //-b
    //================================================ Verbosity
    else if (argument.find("-v") != std::string::npos)
    {
      if ((i + 1) >= argc)
      {
        std::cerr << "Invalid option used with command line argument "
                     "-v. Options are 0,1 or 2."
                  << std::endl;
        chi::Exit(EXIT_FAILURE);
      }
      else
      {
        std::string v_option(argv[i + 1]);
        try
        {
          int level = std::stoi(v_option);
          chi::log.SetVerbosity(level);
        }
        catch (const std::invalid_argument& e)
        {
          std::cerr << "Invalid option used with command line argument "
                       "-v. Options are 0,1 or 2."
                    << std::endl;
          chi::Exit(EXIT_FAILURE);
        }
      }

    } //-v
    else if ((argument.find('=') == std::string::npos) and (!input_file_found))
    {
      chi::run_time::input_file_name_ = argument;
      input_file_found = true;
      chi::run_time::sim_option_interactive_ = false;
    } // no =
    else if (argument.find('=') != std::string::npos)
    {
      chi::console.GetCommandBuffer().push_back(argument);
    } //=

  } // for argument

  if (chi::run_time::dump_registry_)
    ChiObjectMaker::GetInstance().DumpRegister();
}

// ############################################### Initialize ChiTech
/**Initializes all necessary items for ChiTech.
\param argc int    Number of arguments supplied.
\param argv char** Array of strings representing each argument.
 */
int chi::Initialize(int argc, char** argv)
{
  int location_id = 0, number_processes = 1;

  MPI_Init(&argc, &argv);                           /* starts MPI */
  MPI_Comm_rank(MPI_COMM_WORLD, &location_id);      /* get cur process id */
  MPI_Comm_size(MPI_COMM_WORLD, &number_processes); /* get num of processes */

  mpi.SetLocationID(location_id);
  mpi.SetProcessCount(number_processes);

  chi::console.LoadRegisteredLuaItems();
  chi::console.PostMPIInfo(location_id,       number_processes);

  run_time::ParseArguments(argc, argv);

  run_time::InitPetSc(argc, argv);

  return 0;
}

/**Initializes PetSc for use by all entities.*/
int chi::run_time::InitPetSc(int argc, char** argv)
{
  PetscErrorCode ierr;
  PetscMPIInt size;

  PetscOptionsInsertString(nullptr, "-error_output_stderr");
  if (not chi::run_time::allow_petsc_error_handler_)
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
void chi::Finalize()
{
  meshhandler_stack.clear();

  surface_mesh_stack.clear();
  logicvolume_stack.clear();
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
int chi::RunInteractive(int argc, char** argv)
{
  if (not chi::run_time::supress_beg_end_timelog_)
  {
    chi::log.Log() << chi_objects::ChiTimer::GetLocalDateTimeString()
                   << " Running ChiTech in interactive-mode with "
                   << chi::mpi.process_count << " processes.";

    chi::log.Log() << "ChiTech version " << GetVersionStr();
  }

  chi::log.Log() << "ChiTech number of arguments supplied: " << argc - 1;

  chi::log.LogAll();

  chi::console.FlushConsole();

  const auto& input_fname = chi::run_time::input_file_name_;

  if (not input_fname.empty())
  {
    try
    {
      chi::console.ExecuteFile(input_fname, argc, argv);
    }
    catch (const std::exception& excp)
    {
      chi::log.LogAllError() << excp.what();
      // No quitting if file execution fails
    }
  }

  chi::console.RunConsoleLoop();

  if (not chi::run_time::supress_beg_end_timelog_)
  {
    chi::log.Log() << "Final program time " << program_timer.GetTimeString();
    chi::log.Log() << chi_objects::ChiTimer::GetLocalDateTimeString()
                   << " ChiTech finished execution.";
  }

  return 0;
}

// ############################################### Batch interface
/**Runs ChiTech in pure batch mode. Start then finish.*/
int chi::RunBatch(int argc, char** argv)
{
  if (not chi::run_time::supress_beg_end_timelog_)
  {
    chi::log.Log() << chi_objects::ChiTimer::GetLocalDateTimeString()
                   << " Running ChiTech in batch-mode with "
                   << chi::mpi.process_count << " processes.";

    chi::log.Log() << "ChiTech version " << GetVersionStr();
  }

  chi::log.Log() << "ChiTech number of arguments supplied: " << argc - 1;

  if (argc <= 1) chi::log.Log() << chi::run_time::command_line_help_string_;
  chi::console.FlushConsole();

#ifndef NDEBUG
  chi::log.Log() << "Waiting...";
  if (chi::mpi.location_id == 0)
    for (int k = 0; k < 2; ++k)
    {
      usleep(1000000);
      chi::log.Log() << k;
    }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  const auto& input_fname = chi::run_time::input_file_name_;
  int error_code = 0;

  if ((not input_fname.empty()) and (not chi::run_time::termination_posted_))
  {
    try
    {
      error_code = chi::console.ExecuteFile(input_fname, argc, argv);
    }
    catch (const std::exception& excp)
    {
      chi::log.LogAllError() << excp.what();
      chi::Exit(EXIT_FAILURE);
    }
  }

  if (not chi::run_time::supress_beg_end_timelog_)
  {
    chi::log.Log() << "\nFinal program time " << program_timer.GetTimeString();
    chi::log.Log() << chi_objects::ChiTimer::GetLocalDateTimeString()
                   << " ChiTech finished execution of "
                   << chi::run_time::input_file_name_;
  }

  return error_code;
}

// ###################################################################
/** Exits the program appropriately.*/
void chi::Exit(int error_code) { MPI_Abort(MPI_COMM_WORLD, error_code); }

// ###################################################################
/** Gets the ChiTech-version string.*/
std::string chi::GetVersionStr() { return PROJECT_VERSION; }
