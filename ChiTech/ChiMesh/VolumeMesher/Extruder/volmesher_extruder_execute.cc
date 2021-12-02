#include "volmesher_extruder.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/SurfaceMesher/surfacemesher.h"
#include "ChiMesh/Region/chi_region.h"
#include "ChiMesh/Boundary/chi_boundary.h"
#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

#include "ChiConsole/chi_console.h"
extern ChiConsole&   chi_console;

#include <iostream>
#include <vector>

//###################################################################
/**Execution... nough said.*/
void chi_mesh::VolumeMesherExtruder::Execute()
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " VolumeMesherExtruder executed. Memory in use = "
    << chi_console.GetMemoryUsageInMB() << " MB"
    << std::endl;

  //================================================== Get the current handler
  auto mesh_handler = chi_mesh::GetCurrentHandler();
  auto region = mesh_handler->region_stack.back();

  //================================================== Loop over all regions
  chi_log.Log(LOG_0VERBOSE_1)
    << "VolumeMesherExtruder: Processing Region"
    << std::endl;

  //=========================================== Create new continuum
  auto grid = chi_mesh::MeshContinuum::New();
  auto temp_grid = chi_mesh::MeshContinuum::New();
  AddContinuumToRegion(grid, *region);

  //================================== Setup layers
  // populates vertex-layers
  chi_log.Log(LOG_0VERBOSE_1)
    << "VolumeMesherExtruder: Setting up layers" << std::endl;
  SetupLayers();

  //================================== Process templates
  if (template_type == TemplateType::SURFACE_MESH)
  {
    throw std::logic_error("VolumeMesherExtruder: Surfacemesh extrusions"
                           " no longer supported.");
  }
  else if (template_type == TemplateType::UNPARTITIONED_MESH)
  {
    chi_log.Log(LOG_0VERBOSE_1)
      << "VolumeMesherExtruder: Processing unpartitioned mesh"
      << std::endl;

    //================================== Get node_z_incr
    node_z_index_incr = template_unpartitioned_mesh->vertices.size();

    //================================== Create baseline polygons in template
    //                                   continuum
    chi_log.Log(LOG_0VERBOSE_1)
      << "VolumeMesherExtruder: Creating template cells" << std::endl;
    CreatePolygonCells(*template_unpartitioned_mesh, temp_grid);
  }

  chi_log.Log(LOG_0VERBOSE_1)
    << "VolumeMesherExtruder: Creating local nodes" << std::endl;
  CreateLocalNodes(*temp_grid, *grid);

  chi_log.Log(LOG_0VERBOSE_1)
    << "VolumeMesherExtruder: Done creating local nodes" << std::endl;

  //================================== Create extruded item_id
  chi_log.Log(LOG_0)
    << "VolumeMesherExtruder: Extruding cells" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  ExtrudeCells(*temp_grid, *grid);

  size_t total_local_cells = grid->local_cells.size();
  size_t total_global_cells = 0;

  MPI_Allreduce(&total_local_cells,
                &total_global_cells,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,
                MPI_COMM_WORLD);

  chi_log.Log(LOG_0)
    << "VolumeMesherExtruder: Cells extruded = "
    << total_global_cells
    << std::endl;

  //================================== Checking partitioning parameters
  if (options.partition_type != KBA_STYLE_XYZ)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Any partitioning scheme other than KBA_STYLE_XYZ is currently not"
         " supported by VolumeMesherExtruder. No worries. There are plans"
         " to develop this support.";
    exit(EXIT_FAILURE);
  }
  if (!options.mesh_global)
  {
    int p_tot = options.partition_x*options.partition_y*options.partition_z;

    if (chi_mpi.process_count != p_tot)
    {
      chi_log.Log(LOG_ALLERROR)
        << "ERROR: Number of processors available ("
        << chi_mpi.process_count << ") does not match amount of processors "
        << "required by surface mesher partitioning parameters ("
        << p_tot << ").";
      exit(EXIT_FAILURE);
    }
  }//if mesh-global

  chi_log.Log(LOG_ALLVERBOSE_1) << "Building local cell indices";

  //================================== Print info
  chi_log.Log(LOG_ALLVERBOSE_1)
    << "### LOCATION[" << chi_mpi.location_id
    << "] amount of local cells="
    << grid->local_cell_glob_indices.size();


  chi_log.Log(LOG_0)
    << "VolumeMesherExtruder: Number of cells in region = "
    << total_global_cells
    << std::endl;

//  chi_log.Log(LOG_0)
//    << "VolumeMesherExtruder: Number of nodes in region = "
//    << grid->vertices.size()
//    << std::endl;
//  grid->vertices.shrink_to_fit();

  MPI_Barrier(MPI_COMM_WORLD);
}