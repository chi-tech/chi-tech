#include "volmesher_predefined2d.h"
#include <iostream>
#include <vector>
#include "../../MeshHandler/chi_meshhandler.h"
#include "../../Region/chi_region.h"
#include "../../Cell/cell_triangle.h"
#include "../../Boundary/chi_boundary.h"
#include "../../SurfaceMesher/surfacemesher.h"
#include "../../../ChiMPI/chi_mpi.h"
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

#include <ChiTimer/chi_timer.h>
extern ChiTimer chi_program_timer;

void chi_mesh::VolumeMesherPredefined2D::Execute()
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " VolumeMesherPredefined2D executed"
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
    //=========================================== Check for interfaces

    //=========================================== Create new continuum
    //chi_mesh::MeshContinuum* remeshed_surfcont = region->mesh_continua.back();
    chi_mesh::MeshContinuum* vol_continuum = new chi_mesh::MeshContinuum;
    region->volume_mesh_continua.push_back(vol_continuum);

    std::vector<chi_mesh::Boundary*>::iterator bndry;
    //=========================================== Find the first boundary that
    //                                            has a surface mesh and execute
    //                                            the meshing
    bool single_surfacemesh_processed = false;

    for (bndry = region->boundaries.begin();
         bndry != region->boundaries.end();
         bndry++)
    {
      if ((*bndry)->initial_mesh_continuum.surface_mesh!= nullptr)
      {
        //================================== Check if a surface has already
        //                                   been processed
        if (single_surfacemesh_processed)
        {
          std::cerr << "ERROR: Only 1 SurfaceMesh Boundary may be specified ";
          std::cerr << "for VolumeMesherPredefined2D.";
          exit(EXIT_FAILURE);
        }
        else
        {single_surfacemesh_processed = true;}

        //================================== Assign reference continuum, if
        //                                   a mesh operation has been performed
        //                                   after loading the surface was
        //                                   loaded the latest mesh will be used
        chi_mesh::MeshContinuum* ref_continuum =
                    &(*bndry)->initial_mesh_continuum;
        if ((*bndry)->mesh_continua.size()>0)
        {
          ref_continuum = (*bndry)->mesh_continua.back();
        }

        //================================== Create node for each vertex
//        std::vector<chi_mesh::Vertex>::iterator vertex;
//        for (vertex = ref_continuum->surface_mesh->vertices.begin();
//             vertex != ref_continuum->surface_mesh->vertices.end();
//             vertex++)
//        {
//          chi_mesh::Node* node = new chi_mesh::Node;
//          *node = (*vertex.base());
//
//          vol_continuum->nodes.push_back(node);
//        }

        //================================== Create cell for each face
        if (this->options.force_polygons)
        {
          this->CreatePolygonCells(ref_continuum->surface_mesh,
                                   vol_continuum);
        }
        else
        {
          this->CreateTriangleCells(ref_continuum->surface_mesh,
                                    vol_continuum);
        }

        //================================== Connect Boundaries
        std::vector<chi_mesh::Cell*>::iterator cell;
        for (cell = vol_continuum->cells.begin();
             cell != vol_continuum->cells.end();
             cell++)
        {
          (*cell)->FindBoundary2D(region);
        }

        //================================== Check all open item_id have
        //                                   boundaries
        int no_boundary_cells=0;
        for (cell = vol_continuum->cells.begin();
             cell != vol_continuum->cells.end();
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
            << vol_continuum->cells.size()
            << " item_id found with no boundary connection.\n";
          //temp_continuum->ExportCellsToPython("Zerror.py");
        }

        //================================== Checking partitioning parameters
        int p_tot = mesh_handler->surface_mesher->partitioning_x*
                    mesh_handler->surface_mesher->partitioning_y;
        if (chi_mpi.process_count != p_tot)
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

        //================================== InitializeAlphaElements local cell indices
        int num_glob_cells=vol_continuum->cells.size();
        for (int c=0; c<num_glob_cells; c++)
        {
          vol_continuum->glob_cell_local_indices.push_back(-1);
          if (vol_continuum->cells[c]->partition_id ==
              chi_mpi.location_id)
          {
            vol_continuum->local_cell_glob_indices.push_back(c);
            int local_cell_index = vol_continuum->local_cell_glob_indices.size()-1;
            vol_continuum->glob_cell_local_indices[c]=local_cell_index;
            vol_continuum->cells[c]->cell_local_id = local_cell_index;
          }
        }
        chi_log.Log(LOG_ALLVERBOSE_1)
        << "### LOCATION[" << chi_mpi.location_id
        << "] amount of local cells="
        << vol_continuum->local_cell_glob_indices.size();


        chi_log.Log(LOG_0)
          << "VolumeMesherPredefined2D["
          << chi_mpi.location_id
          << "]: Number of cells in region = "
          << vol_continuum->cells.size()
          << std::endl;

        chi_log.Log(LOG_0)
          << "VolumeMesherPredefined2D["
          << chi_mpi.location_id
          << "]: Number of nodes in region = "
          << vol_continuum->nodes.size()
          << std::endl;


//        if (chi_mpi.location_id == 0)
//        {
//          vol_continuum->ExportCellsToPython("SurfaceMesh.py");
//        }
      } //if surface mesh
    } //for boundaries
  } //for regions

  MPI_Barrier(MPI_COMM_WORLD);
}