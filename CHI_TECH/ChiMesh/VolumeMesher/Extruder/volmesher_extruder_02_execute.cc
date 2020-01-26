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
  std::vector<chi_mesh::Region*>::iterator region_iter;
  for (region_iter = mesh_handler->region_stack.begin();
       region_iter != mesh_handler->region_stack.end();
       region_iter++)
  {
    chi_mesh::Region* region = *region_iter;

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
    region->volume_mesh_continua.push_back(temp_grid);
    region->volume_mesh_continua.push_back(grid);

    //=========================================== Perform the operation
    bool single_surfacemesh_processed = false;

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
        bool delete_surface_mesh_elements = true;
        CreatePolygonCells(ref_continuum->surface_mesh, temp_grid,
                           delete_surface_mesh_elements);
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

        chi_log.Log(LOG_0)
          << "VolumeMesherExtruder: Cells extruded = "
          << grid->cells.size()
          << std::endl;



        //================================== Clean-up temporary continuum
        for (auto vert : temp_grid->nodes) delete vert;
        for (auto pcell : temp_grid->cells) delete pcell;
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
        int num_glob_cells=grid->cells.size();
        for (int c=0; c<num_glob_cells; c++)
        {
          grid->glob_cell_local_indices.push_back(-1);

          if (grid->cells[c] == nullptr)
            continue;

          if ((grid->cells[c]->partition_id == chi_mpi.location_id) ||
              (options.mesh_global))
          {
            grid->local_cell_glob_indices.push_back(c);
            int local_cell_index = grid->local_cell_glob_indices.size() - 1;
            grid->glob_cell_local_indices[c]=local_cell_index;

            grid->cells[c]->cell_local_id = local_cell_index;
          }
        }

        chi_log.Log(LOG_ALLVERBOSE_1)
          << "### LOCATION[" << chi_mpi.location_id
          << "] amount of local cells="
          << grid->local_cell_glob_indices.size();


        chi_log.Log(LOG_0)
          << "VolumeMesherExtruder: Number of cells in region = "
          << grid->cells.size()
          << std::endl;
        grid->cells.shrink_to_fit();

        chi_log.Log(LOG_0)
          << "VolumeMesherExtruder: Number of nodes in region = "
          << grid->nodes.size()
          << std::endl;
        grid->nodes.shrink_to_fit();


      }//if surface mesh
    }//for bndry
  }//for regions

  MPI_Barrier(MPI_COMM_WORLD);
}