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


//###################################################################
/**Execution... nough said.*/
void chi_mesh::VolumeMesherExtruder::Execute()
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " VolumeMesherExtruder executed"
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
    chi_mesh::Boundary* bot_boundary = new chi_mesh::Boundary;
    chi_mesh::Boundary* top_boundary = new chi_mesh::Boundary;
    region->boundaries.push_back(bot_boundary);
    region->boundaries.push_back(top_boundary);
    bot_boundary_index = region->boundaries.size()-2;
    top_boundary_index = region->boundaries.size()-1;


    //=========================================== Check for interfaces

    //=========================================== Create new continuum
    //chi_mesh::MeshContinuum* remeshed_surfcont = region->mesh_continua.back();
    auto vol_continuum = new chi_mesh::MeshContinuum;
    auto temp_continuum = new chi_mesh::MeshContinuum;
    region->volume_mesh_continua.push_back(temp_continuum);
    region->volume_mesh_continua.push_back(vol_continuum);

    std::vector<chi_mesh::Boundary*>::iterator bndry;
    //=========================================== Perform the operation
    bool single_surfacemesh_processed = false;

    for (bndry = region->boundaries.begin();
         bndry != region->boundaries.end();
         bndry++)
    {
      if ((*bndry)->initial_mesh_continuum.surface_mesh!= nullptr)
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
        chi_mesh::MeshContinuum* ref_continuum = &(*bndry)->initial_mesh_continuum;
        if ((*bndry)->mesh_continua.size()>0)
        {
          ref_continuum = (*bndry)->mesh_continua.back();
        }

        //We now have the surface we want to extrude
        chi_log.Log(LOG_0VERBOSE_1)
          << "VolumeMesherExtruder: Setting up layers"
          << std::endl;
        SetupLayers();


        //================================== Create nodes
        chi_log.Log(LOG_0VERBOSE_1)
          << "VolumeMesherExtruder: Creating nodes"
          << std::endl;

        node_z_index_incr = 0;
        for (int iz=0; iz<vertex_layers.size(); iz++)
        {
          std::vector<chi_mesh::Vertex>::iterator vertex;
          for (vertex = ref_continuum->surface_mesh->vertices.begin();
               vertex != ref_continuum->surface_mesh->vertices.end();
               vertex++)
          {
            chi_mesh::Node* node = new chi_mesh::Node;
            *node = (*vertex.base());

            node->z = vertex_layers[iz];

            vol_continuum->nodes.push_back(node);
          }
          if (iz==0)
          {
            node_z_index_incr = vol_continuum->nodes.size();
          }
        }

        //================================== Create baseline polygons in template
        //                                   continuum
        chi_log.Log(LOG_0VERBOSE_1)
          << "VolumeMesherExtruder: Creating template cells"
          << std::endl;
        CreatePolygonCells(ref_continuum->surface_mesh, temp_continuum);


        //================================== Connect template Boundaries
        chi_log.Log(LOG_0VERBOSE_1)
          << "VolumeMesherExtruder: Connecting boundaries"
          << std::endl;
        std::vector<chi_mesh::Cell*>::iterator cell;
        for (cell = temp_continuum->cells.begin();
             cell != temp_continuum->cells.end();
             cell++)
        {
          (*cell)->FindBoundary2D(region);
        }

        //================================== Check all open item_id of template
        //                                   have boundaries
        int no_boundary_cells=0;
        for (cell = temp_continuum->cells.begin();
             cell != temp_continuum->cells.end();
             cell++)
        {
          if (!(*cell)->CheckBoundary2D())
          {
            no_boundary_cells++;
          }

        }
        if (no_boundary_cells>0)
        {
          chi_log.Log(LOG_ALLVERBOSE_1)
          << "A total of "
          << no_boundary_cells
          << " out of "
          << temp_continuum->cells.size()
          << " item_id found with no boundary connection.\n";
          //temp_continuum->ExportCellsToPython("Zerror.py");
        }
        //================================== Create extruded item_id
        chi_log.Log(LOG_0)
          << "VolumeMesherExtruder: Extruding cells" << std::endl;

        ExtrudeCells(temp_continuum, vol_continuum);

        chi_log.Log(LOG_0)
          << "VolumeMesherExtruder: Cells extruded = "
          << vol_continuum->cells.size()
          << std::endl;


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
        }


        //================================== Initialize local cell indices
        int num_glob_cells=vol_continuum->cells.size();
        for (int c=0; c<num_glob_cells; c++)
        {
          vol_continuum->glob_cell_local_indices.push_back(-1);
          if ((vol_continuum->cells[c]->partition_id == chi_mpi.location_id) ||
              (options.mesh_global))
          {
            vol_continuum->local_cell_glob_indices.push_back(c);
            int local_cell_index =
              vol_continuum->local_cell_glob_indices.size()-1;
            vol_continuum->glob_cell_local_indices[c]=local_cell_index;

            vol_continuum->cells[c]->cell_local_id = local_cell_index;
          }
        }
        chi_log.Log(LOG_ALLVERBOSE_1)
          << "### LOCATION[" << chi_mpi.location_id
          << "] amount of local cells="
          << vol_continuum->local_cell_glob_indices.size();


        chi_log.Log(LOG_0)
          << "VolumeMesherExtruder: Number of cells in region = "
          << vol_continuum->cells.size()
          << std::endl;
        vol_continuum->cells.shrink_to_fit();

        chi_log.Log(LOG_0)
          << "VolumeMesherExtruder: Number of nodes in region = "
          << vol_continuum->nodes.size()
          << std::endl;
        vol_continuum->nodes.shrink_to_fit();


      }//if surface mesh
    }//for bndry
  }//for regions

  MPI_Barrier(MPI_COMM_WORLD);
}