#include "volmesher_linemesh1d.h"

#include "../../MeshHandler/chi_meshhandler.h"
#include "../../Region/chi_region.h"
#include "../../Boundary/chi_boundary.h"
#include "../../Cell/cell_slab.h"
#include "../../SurfaceMesher/surfacemesher.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI chi_mpi;
extern ChiLog chi_log;

#include <ChiTimer/chi_timer.h>
extern ChiTimer chi_program_timer;

//###################################################################
/**Executes the One dimensional mesher process.*/
void chi_mesh::VolumeMesherLinemesh1D::Execute()
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " VolumeMesherLinemesh1D executed"
    << std::endl;

  //================================================== Get the current handler
  chi_mesh::MeshHandler* mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::SurfaceMesher* surf_mesher = mesh_handler->surface_mesher;

  //================================================== Loop over all regions
  std::vector<chi_mesh::Region*>::iterator region_iter;
  for (region_iter = mesh_handler->region_stack.begin();
       region_iter != mesh_handler->region_stack.end();
       region_iter++)
  {
    chi_mesh::Region* region = *region_iter;

    chi_log.Log(LOG_0VERBOSE_1)
      << "VolumeMesherLinemesh1D: Processing Region"
      << std::endl;

    //=========================================== Create new continuum
    chi_mesh::MeshContinuum* vol_continuum = new chi_mesh::MeshContinuum;
    region->volume_mesh_continua.push_back(vol_continuum);

    std::vector<chi_mesh::Boundary*>::iterator bndry;
    //=========================================== Perform the operation
    bool single_linemesh_processed = false;

    for (bndry = region->boundaries.begin();
         bndry != region->boundaries.end();
         bndry++)
    {
      if ((*bndry)->initial_mesh_continuum.line_mesh != nullptr)
      {
        //================================== Check for duplicate linemesh
        if (single_linemesh_processed)
        {
          std::cerr << "ERROR: Only 1 LineMesh Boundary may be specified ";
          std::cerr << "for VolumeMesherLinemesh1D.";
          exit(EXIT_FAILURE);
        }
        else
        {single_linemesh_processed = true;}

        chi_mesh::LineMesh* line_mesh =
          (*bndry)->initial_mesh_continuum.line_mesh;

        //================================== Populate nodes
        for (int v=0; v<line_mesh->vertices.size(); v++)
        {
          vol_continuum->nodes.push_back(&line_mesh->vertices[v]);
        }
        num_slab_cells = line_mesh->vertices.size()-1;

        //================================== Create cells from line mesh
        int cell_count = -1;
        for (int v=0; v<(line_mesh->vertices.size()-1); v++)
        {
          cell_count++;
          chi_mesh::CellSlab* slab = new chi_mesh::CellSlab;
          slab->cell_global_id = vol_continuum->cells.size();

          //====================== Populate basic data
          slab->v_indices[0] = v;
          slab->v_indices[1] = v+1;

          slab->edges[0] = cell_count-1;
          slab->edges[1] = cell_count+1;

          if (v == (line_mesh->vertices.size()-2))
            slab->edges[1] = -2;

          //====================== Compute centroid
          chi_mesh::Vertex v0 = line_mesh->vertices[v];
          chi_mesh::Vertex v1 = line_mesh->vertices[v+1];

          slab->centroid = (v0+v1)/2.0;

          //====================== Compute normals
          chi_mesh::Vector n = (v1-v0)/(v1-v0).Norm();
          slab->face_normals[0] = chi_mesh::Vector(0.0,0.0,-1.0);
          slab->face_normals[1] = chi_mesh::Vector(0.0,0.0, 1.0);

          slab->xyz_partition_indices = GetCellXYZPartitionID(slab);

          int xi,yi,zi;
          xi = std::get<0>(slab->xyz_partition_indices);
          yi = std::get<1>(slab->xyz_partition_indices);
          zi = std::get<2>(slab->xyz_partition_indices);
          int px,py;
          px = surf_mesher->partitioning_x;
          py = surf_mesher->partitioning_y;
          slab->partition_id = zi*px*py + yi*px + xi;

          vol_continuum->cells.push_back(slab);

        }//for interval

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
          << "VolumeMesherLinemesh1D: Number of cells in region = "
          << vol_continuum->cells.size()
          << std::endl;
        vol_continuum->cells.shrink_to_fit();

        chi_log.Log(LOG_0)
          << "VolumeMesherLinemesh1D: Number of nodes in region = "
          << vol_continuum->nodes.size()
          << std::endl;
        vol_continuum->nodes.shrink_to_fit();

      }//if linemesh
    }//for bndry
  }//for regions

  MPI_Barrier(MPI_COMM_WORLD);
}