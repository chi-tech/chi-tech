#include "volmesher_predefined2d.h"
#include <iostream>
#include <vector>
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/Region/chi_region.h"
#include "ChiMesh/Boundary/chi_boundary.h"
#include "ChiMesh/SurfaceMesher/surfacemesher.h"
#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

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
  int total_global_cells = 0;
  for (auto region : mesh_handler->region_stack)
  {
    //=========================================== Create new continuum
    auto grid = new chi_mesh::MeshContinuum;
    AddContinuumToRegion(grid, *region);

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

        //================================== Create cell for each face
        auto temp_grid = new chi_mesh::MeshContinuum;
        this->CreatePolygonCells(ref_continuum->surface_mesh, temp_grid);
        GridFilterGhosts(temp_grid,grid);

//        //================================== Connect Boundaries
//        for (auto& cell : grid->local_cells)
//          cell.FindBoundary2D(region);

//        //================================== Check all open item_id have
//        //                                   boundaries
//        int no_boundary_cells=0;
//        for (auto cell : grid->local_cells)
//          if (!cell.CheckBoundary2D())
//            no_boundary_cells++;

        int total_local_cells = grid->local_cells.size();

        MPI_Allreduce(&total_local_cells,
                      &total_global_cells,
                      1,
                      MPI_INT,
                      MPI_SUM,
                      MPI_COMM_WORLD);


//        if (no_boundary_cells>0)
//        {
//          chi_log.Log(LOG_ALLVERBOSE_1)
//            << "A total of "
//            << no_boundary_cells
//            << " out of "
//            << total_global_cells
//            << " item_id found with no boundary connection.\n";
//          //temp_continuum->ExportCellsToPython("Zerror.py");
//        }

        //================================== Checking partitioning parameters
        if (!options.mesh_global)
        {
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
        }


        //================================== InitializeAlphaElements local cell indices
        chi_log.Log(LOG_ALLVERBOSE_1)
          << "### LOCATION[" << chi_mpi.location_id
          << "] amount of local cells="
          << grid->local_cell_glob_indices.size();


        chi_log.Log(LOG_0)
          << "VolumeMesherPredefined2D["
          << chi_mpi.location_id
          << "]: Number of cells in region = "
          << total_global_cells
          << std::endl;

        chi_log.Log(LOG_0)
          << "VolumeMesherPredefined2D["
          << chi_mpi.location_id
          << "]: Number of nodes in region = "
          << grid->vertices.size()
          << std::endl;

      } //if surface mesh
    } //for boundaries
  } //for regions

  MPI_Barrier(MPI_COMM_WORLD);
}