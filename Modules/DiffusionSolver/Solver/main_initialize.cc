#include "diffusion_solver.h"



#include <ChiTimer/chi_timer.h>
#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;
extern ChiTimer chi_program_timer;

//###################################################################
/**Initializes the diffusion solver using the PETSc library.*/
int chi_diffusion::Solver::Initialize(bool verbose)
{
  chi_log.Log(LOG_0) << "\n"
                     << chi_program_timer.GetTimeString() << " "
                     << solver_name << ": Initializing Diffusion solver PETSc";
  this->verbose_info = verbose;

  if (regions.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_diffusion::Solver::Initialize: No region added to solver.";
    exit(EXIT_FAILURE);
  }

  chi_mesh::Region* region = regions.back();
  grid = region->GetGrid();



  if (not common_items_initialized)
    InitializeCommonItems(); //Mostly boundaries

  ChiTimer t_init; t_init.Reset();

  //================================================== Initialize discretization
  //                                                   method
  if (fem_method == PWLC)
    InitializePWLC(verbose);
  else if (fem_method == PWLD_MIP)
    InitializePWLD(verbose);
  else if (fem_method == PWLD_MIP_GRPS)
    InitializePWLDGroups(verbose);
  else if (fem_method == PWLD_MIP_GAGG)
    InitializePWLDGrpAgg(verbose);
  else
  {
    chi_log.Log(LOG_0)
      << "Diffusion Solver: Finite Element Discretization "
         "method not specified.";
    exit(EXIT_FAILURE);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString() << " "
    << solver_name << ": Diffusion Solver initialization time "
    << t_init.GetTime()/1000.0 << std::endl;

  return false;
}