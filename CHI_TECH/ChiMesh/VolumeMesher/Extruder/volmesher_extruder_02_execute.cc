#include "volmesher_extruder.h"
#include <iostream>
#include <vector>
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/SurfaceMesher/surfacemesher.h"
#include "ChiMesh/Region/chi_region.h"
#include "ChiMesh/Boundary/chi_boundary.h"
#include <ChiMPI/chi_mpi.h>

#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

#include <ChiTimer/chi_timer.h>
extern ChiTimer chi_program_timer;

#include <ChiConsole/chi_console.h>
extern ChiConsole  chi_console;


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
  chi_mesh::MeshHandler* mesh_handler = chi_mesh::GetCurrentHandler();

  //================================================== Loop over all regions
  bool single_surfacemesh_processed = false;
  int total_global_cells = 0;
  for (auto region : mesh_handler->region_stack)
  {
    chi_log.Log(LOG_0VERBOSE_1)
      << "VolumeMesherExtruder: Processing Region"
      << std::endl;

    //=========================================== Add top and bottom boundaries
    auto bot_boundary = new chi_mesh::Boundary;
    auto top_boundary = new chi_mesh::Boundary;
    region->boundaries.push_back(bot_boundary);
    region->boundaries.push_back(top_boundary);
    bot_boundary_index = region->boundaries.size()-2;
    top_boundary_index = region->boundaries.size()-1;

    //=========================================== Check for interfaces

    //=========================================== Create new continuum
    auto grid = new chi_mesh::MeshContinuum;
    auto temp_grid = new chi_mesh::MeshContinuum;
    AddContinuumToRegion(grid, *region);

    //=========================================== Look over boundaries
    for (auto bndry : region->boundaries)
    {
      if (bndry->initial_mesh_continuum.surface_mesh!= nullptr)
      {
        chi_log.Log(LOG_0VERBOSE_1)
          << "VolumeMesherExtruder: Processing surface mesh"
          << std::endl;

        //================================== Check for duplicate surface
        if (single_surfacemesh_processed)
        {
          std::cerr << "ERROR: Only 1 SurfaceMesh Boundary may be specified ";
          std::cerr << "for VolumeMesherExtruder.";
          exit(EXIT_FAILURE);
        }
        else
        {single_surfacemesh_processed = true;}

        //================================== Assign reference continuum
        chi_mesh::MeshContinuum* ref_continuum = &bndry->initial_mesh_continuum;
        if (not bndry->mesh_continua.empty())
          ref_continuum = bndry->mesh_continua.back();

        //We now have the surface we want to extrude
        chi_log.Log(LOG_0VERBOSE_1)
          << "VolumeMesherExtruder: Setting up layers"
          << std::endl;
        SetupLayers();

        //================================== Get node_z_incr
        node_z_index_incr = ref_continuum->surface_mesh->vertices.size();

        //================================== Create baseline polygons in template
        //                                   continuum
        chi_log.Log(LOG_0VERBOSE_1)
          << "VolumeMesherExtruder: Creating template cells"
          << std::endl;
        const bool DELETE_SURFACE_MESH_ELEMENTS = true;
        const bool FORCE_LOCAL = true;
        CreatePolygonCells(ref_continuum->surface_mesh,
                           temp_grid,
                           DELETE_SURFACE_MESH_ELEMENTS,
                           FORCE_LOCAL);
        delete ref_continuum->surface_mesh;

        chi_log.Log(LOG_0VERBOSE_1)
          << "VolumeMesherExtruder: Creating local nodes"
          << std::endl;
        CreateLocalAndBoundaryNodes(temp_grid,grid);

        chi_log.Log(LOG_0VERBOSE_1)
          << "VolumeMesherExtruder: Done creating local nodes"
          << std::endl;

        //================================== Create extruded item_id
        chi_log.Log(LOG_0)
          << "VolumeMesherExtruder: Extruding cells" << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        ExtrudeCells(temp_grid, grid);

        int total_local_cells = grid->local_cells.size();

        MPI_Allreduce(&total_local_cells,
                      &total_global_cells,
                      1,
                      MPI_INT,
                      MPI_SUM,
                      MPI_COMM_WORLD);

        chi_log.Log(LOG_0)
          << "VolumeMesherExtruder: Cells extruded = "
          << total_global_cells
          << std::endl;



        //================================== Clean-up temporary continuum
        for (auto vert : temp_grid->vertices) delete vert;
        delete temp_grid;

        //================================== Checking partitioning parameters
        if (!options.mesh_global)
        {
          int p_tot = mesh_handler->surface_mesher->partitioning_x*
                      mesh_handler->surface_mesher->partitioning_y*
                      options.partition_z;
          chi_log.Log(LOG_ALLVERBOSE_2)
            << "Processes called for " << p_tot
            << ". Processes supplied " << chi_mpi.process_count;

          if ((chi_mpi.process_count != p_tot) /*&& (p_tot != 0)*/)
          {
            chi_log.Log(LOG_ALLERROR) <<
                                      "ERROR: Number of processors available ("
                                      << chi_mpi.process_count <<
                                      ") does not match amount of processors "
                                      "required by surface"
                                      " mesher partitioning parameters ("
                                      << p_tot <<
                                      ").";
            exit(EXIT_FAILURE);
          }
        }//if mesh-global

        chi_log.Log(LOG_ALLVERBOSE_1) << "Building local cell indices";

        //================================== Initialize local cell indices
        chi_log.Log(LOG_ALLVERBOSE_1)
          << "### LOCATION[" << chi_mpi.location_id
          << "] amount of local cells="
          << grid->local_cell_glob_indices.size();


        chi_log.Log(LOG_0)
          << "VolumeMesherExtruder: Number of cells in region = "
          << total_global_cells
          << std::endl;

        chi_log.Log(LOG_0)
          << "VolumeMesherExtruder: Number of nodes in region = "
          << grid->vertices.size()
          << std::endl;
        grid->vertices.shrink_to_fit();


      }//if surface mesh
    }//for bndry
  }//for regions

  if (not single_surfacemesh_processed)
  {
    chi_log.Log(LOG_ALLERROR)
      << "VolumeMesherExtruder: No surface mesh was processed for any region."
         " Use \"chiRegionAddSurfaceBoundary\" to add a surface to the region.";
    exit(EXIT_FAILURE);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}