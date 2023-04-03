#include "volmesher_extruder.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/SurfaceMesher/surfacemesher.h"
#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"


#include "ChiTimer/chi_timer.h"


#include "ChiConsole/chi_console.h"

#include <iostream>

//###################################################################
/**Execution... nough said.*/
void chi_mesh::VolumeMesherExtruder::Execute()
{
  chi::log.Log()
    << chi::program_timer.GetTimeString()
    << " VolumeMesherExtruder executed. Memory in use = "
    << chi_objects::ChiConsole::GetMemoryUsageInMB() << " MB"
    << std::endl;

  //================================================== Loop over all regions
  chi::log.Log0Verbose1()
    << "VolumeMesherExtruder: Processing Region"
    << std::endl;

  //=========================================== Create new continuum
  auto grid = chi_mesh::MeshContinuum::New();
  auto temp_grid = chi_mesh::MeshContinuum::New();

  SetContinuum(grid);
  SetGridAttributes(DIMENSION_3 | EXTRUDED);

  //================================== Setup layers
  // populates vertex-layers
  chi::log.Log0Verbose1()
    << "VolumeMesherExtruder: Setting up layers" << std::endl;
  SetupLayers();

  //================================== Process templates
  if (template_type_ == TemplateType::UNPARTITIONED_MESH)
  {
    chi::log.Log0Verbose1()
      << "VolumeMesherExtruder: Processing unpartitioned mesh"
      << std::endl;

    //================================== Get node_z_incr
    node_z_index_incr_ = template_unpartitioned_mesh_->GetVertices().size();

    //================================== Create baseline polygons in template
    //                                   continuum
    chi::log.Log0Verbose1()
      << "VolumeMesherExtruder: Creating template cells" << std::endl;
    CreatePolygonCells(*template_unpartitioned_mesh_, temp_grid);

    grid->GetBoundaryIDMap() =
      template_unpartitioned_mesh_->GetMeshOptions().boundary_id_map;
    temp_grid->GetBoundaryIDMap() =
      template_unpartitioned_mesh_->GetMeshOptions().boundary_id_map;
  }

  chi::log.Log0Verbose1()
    << "VolumeMesherExtruder: Creating local nodes" << std::endl;
  CreateLocalNodes(*temp_grid, *grid);

  chi::log.Log0Verbose1()
    << "VolumeMesherExtruder: Done creating local nodes" << std::endl;

  //================================== Insert top and bottom boundary
  //                                   id map
  auto& grid_bndry_id_map = grid->GetBoundaryIDMap();
  zmax_bndry_id = grid->MakeBoundaryID("ZMAX");
  grid_bndry_id_map[zmax_bndry_id] = "ZMAX";
  zmin_bndry_id = grid->MakeBoundaryID("ZMIN");
  grid_bndry_id_map[zmin_bndry_id] = "ZMIN";


  //================================== Create extruded item_id
  chi::log.Log()
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

  chi::log.Log()
    << "VolumeMesherExtruder: Cells extruded = "
    << total_global_cells
    << std::endl;

  //================================== Checking partitioning parameters
  if (options.partition_type != KBA_STYLE_XYZ)
  {
    chi::log.LogAllError()
      << "Any partitioning scheme other than KBA_STYLE_XYZ is currently not"
         " supported by VolumeMesherExtruder. No worries. There are plans"
         " to develop this support.";
   chi::Exit(EXIT_FAILURE);
  }
  if (!options.mesh_global)
  {
    int p_tot = options.partition_x*options.partition_y*options.partition_z;

    if (chi::mpi.process_count != p_tot)
    {
      chi::log.LogAllError()
        << "ERROR: Number of processors available ("
        << chi::mpi.process_count << ") does not match amount of processors "
        << "required by surface mesher partitioning parameters ("
        << p_tot << ").";
     chi::Exit(EXIT_FAILURE);
    }
  }//if mesh-global

  chi::log.LogAllVerbose1() << "Building local cell indices";

  //================================== Print info
  chi::log.LogAllVerbose1()
    << "### LOCATION[" << chi::mpi.location_id
    << "] amount of local cells="
    << grid->local_cells.size();


  chi::log.Log()
    << "VolumeMesherExtruder: Number of cells in region = "
    << total_global_cells
    << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
}