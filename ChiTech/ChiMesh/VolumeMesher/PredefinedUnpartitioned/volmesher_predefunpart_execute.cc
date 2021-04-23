#include "volmesher_predefunpart.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/SurfaceMesher/surfacemesher.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Region//chi_region.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

#include <ChiTimer/chi_timer.h>
extern ChiTimer chi_program_timer;

#include <ChiConsole/chi_console.h>
extern ChiConsole&   chi_console;

//###################################################################
/**Executes the predefined3D mesher.*/
void chi_mesh::VolumeMesherPredefinedUnpartitioned::Execute()
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " VolumeMesherPredefinedUnpartitioned executing. Memory in use = "
    << chi_console.GetMemoryUsageInMB() << " MB"
    << std::endl;

  //======================================== Get the current handler
  auto mesh_handler = chi_mesh::GetCurrentHandler();

  //======================================== Check empty region list
  if (mesh_handler->region_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "VolumeMesherPredefinedUnpartitioned: No region added.";
    exit(EXIT_FAILURE);
  }

  //======================================== Check unpartitioned mesh available
  if (mesh_handler->unpartitionedmesh_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "VolumeMesherPredefinedUnpartitioned: "
         "No unpartitioned mesh to operate on.";
    exit(EXIT_FAILURE);
  }

  //======================================== Check partitioning params
  if (options.partition_type == KBA_STYLE_XYZ)
  {
    int Px = mesh_handler->volume_mesher->options.partition_x;
    int Py = mesh_handler->volume_mesher->options.partition_y;
    int Pz = mesh_handler->volume_mesher->options.partition_z;

    int desired_process_count = Px*Py*Pz;

    if (desired_process_count != chi_mpi.process_count)
    {
      chi_log.Log(LOG_ALLERROR)
        << "ERROR: Number of processors available ("
        << chi_mpi.process_count <<
        ") does not match amount of processors "
        "required by partitioning parameters ("
        << desired_process_count << ").";
      exit(EXIT_FAILURE);
    }
  }

  //======================================== Get unpartitioned mesh
  auto umesh = mesh_handler->unpartitionedmesh_stack.back();

  chi_log.Log(LOG_0) << "Computed centroids";
  MPI_Barrier(MPI_COMM_WORLD);


  //======================================== Apply partitioning scheme
  auto grid = chi_mesh::MeshContinuum::New();

  if (options.partition_type == PartitionType::KBA_STYLE_XYZ)
    KBA(umesh, grid);
  else
    PARMETIS(umesh,grid);

  chi_log.Log(LOG_0) << "Cells loaded.";
  MPI_Barrier(MPI_COMM_WORLD);

  AddContinuumToRegion(grid, *mesh_handler->region_stack.back());


  //======================================== Concluding messages
  chi_log.Log(LOG_0)
    << "VolumeMesherPredefinedUnpartitioned: Number of nodes in region = "
    << grid->vertices.size()
    << std::endl;
  grid->vertices.shrink_to_fit();

  chi_log.Log(LOG_ALLVERBOSE_1)
    << "### LOCATION[" << chi_mpi.location_id
    << "] amount of local cells="
    << grid->local_cell_glob_indices.size();

  uint64_t total_local_cells = grid->local_cells.size();
  uint64_t total_global_cells = 0;

  MPI_Allreduce(&total_local_cells,
                &total_global_cells,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD);

  chi_log.Log(LOG_0)
    << "VolumeMesherPredefinedUnpartitioned: Cells created = "
    << total_global_cells
    << std::endl;

}