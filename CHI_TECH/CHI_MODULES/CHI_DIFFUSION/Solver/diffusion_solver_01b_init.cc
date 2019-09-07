#include "diffusion_solver.h"

#include"../../../CHI_MESH/CHI_CELL/cell_slab.h"
#include"../../../CHI_MESH/CHI_CELL/cell_triangle.h"
#include "../../../CHI_MESH/CHI_CELL/cell_polygon.h"
#include "../../../CHI_MESH/CHI_CELL/cell_polyhedron.h"
#include "../../../CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../../CHI_MESH/CHI_VOLUMEMESHER/chi_volumemesher.h"

#include <ChiTimer/chi_timer.h>
#include <chi_mpi.h>
#include <chi_log.h>

extern CHI_MPI chi_mpi;
extern CHI_LOG chi_log;
extern ChiTimer chi_program_timer;

#include<fstream>
#include <unistd.h>

//###################################################################
/**Initializes the diffusion solver using the PETSc library.*/
int chi_diffusion::Solver::Initialize(bool verbose)
{
  chi_log.Log(LOG_0) << "\n"
                     << chi_program_timer.GetTimeString() << " "
                     << solver_name << ": Initializing Diffusion solver PETSc";
  this->verbose = verbose;

  if (not common_items_initialized)
    InitializeCommonItems();

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

  //================================================== Compute cell matrices
  chi_mesh::Region*              aregion = this->regions.back();
  chi_mesh::MeshContinuum* vol_continuum =
    aregion->volume_mesh_continua.back();

  if (verbose)
    chi_log.Log(LOG_0) << "Computing cell matrices";
  this->discretization->AddViewOfLocalContinuum(
    vol_continuum,
    vol_continuum->local_cell_glob_indices.size(),
    vol_continuum->local_cell_glob_indices.data());
  MPI_Barrier(MPI_COMM_WORLD);




  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString() << " "
    << solver_name << ": Diffusion Solver initialization time "
    << t_init.GetTime()/1000.0 << std::endl;

  return false;
}