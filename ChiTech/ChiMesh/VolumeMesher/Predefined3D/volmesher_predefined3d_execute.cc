#include "volmesher_predefined3d.h"

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
void chi_mesh::VolumeMesherPredefined3D::Execute()
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " VolumeMesherPredefined3D executing. Memory in use = "
    << chi_console.GetMemoryUsageInMB() << " MB"
    << std::endl;

  //================================================== Get the current handler
  auto mesh_handler = chi_mesh::GetCurrentHandler();

  //======================================== Check empty region list
  if (mesh_handler->region_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "VolumeMesherPredefined3D: No region added.";
    exit(EXIT_FAILURE);
  }

  //======================================== Check empty unpartitioned
  //                                                   list
  if (mesh_handler->unpartitionedmesh_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "VolumeMesherPredefined3D: No unpartitioned mesh to operate on.";
    exit(EXIT_FAILURE);
  }

  //======================================== Check paritioning params
  int Px = 1;
  int Py = 1;
  int Pz = 1;
  if (options.partition_type == KBA_STYLE_XYZ)
  {
    Px = mesh_handler->surface_mesher->partitioning_x;
    Py = mesh_handler->surface_mesher->partitioning_y;
    Pz = mesh_handler->volume_mesher->options.partition_z;

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

  //======================================== Populate vertex
  //                                                   subscriptionns
  auto umesh = mesh_handler->unpartitionedmesh_stack.back();
  std::vector<std::set<int>> vertex_subs(umesh->vertices.size());
  int c=-1;
  for (auto cell : umesh->raw_cells)
  {
    ++c;
    for (auto vid : cell->vertex_ids)
      vertex_subs[vid].insert(c);
  }

  int num_bndry_faces = 0;
  for (auto cell : umesh->raw_cells)
    for (auto& face : cell->faces)
      if (face.neighbor < 0) ++num_bndry_faces;

  chi_log.Log(LOG_0) << "Number of bndry faces: " << num_bndry_faces;

  //======================================== Establish connectivity
  std::vector<std::set<int>> cells_to_search;
  c=-1;
  for (auto cell : umesh->raw_cells)
  {
    ++c;
    for (auto& face : cell->faces)
    {
      if (face.neighbor >= 0) continue;

      cells_to_search.clear();
      for (auto cfvid : face.vertex_ids)
        cells_to_search.push_back(vertex_subs[cfvid]);

      std::set<int> cfvids(face.vertex_ids.begin(),
                           face.vertex_ids.end());

      bool stop_searching = false;
      for (auto& adj_cell_id_set : cells_to_search)
      {
        for (auto& adj_cell_id : adj_cell_id_set)
        {
          if (adj_cell_id == c) continue;
          auto adj_cell = umesh->raw_cells[adj_cell_id];

          bool adj_cell_matches = false;
          for (auto& aface : adj_cell->faces)
          {
            std::set<int> afvids(aface.vertex_ids.begin(),
                                 aface.vertex_ids.end());

            if (afvids == cfvids)
            {
              adj_cell_matches = true;
              break;
            }
          }

          if (adj_cell_matches)
          {
            face.neighbor = adj_cell_id;
            stop_searching = true;
          }

          if (stop_searching) break;
        }//cell id
        if (stop_searching) break;
      }//cell id set
    }//for face
  }//for cell

  num_bndry_faces = 0;
  for (auto cell : umesh->raw_cells)
    for (auto& face : cell->faces)
      if (face.neighbor < 0) ++num_bndry_faces;

  chi_log.Log(LOG_0) << "Number of bndry faces: " << num_bndry_faces;

  //======================================== Compute centroid
  for (auto cell : umesh->raw_cells)
  {
    cell->centroid = chi_mesh::Vertex(0.0,0.0,0.0);
    for (auto vid : cell->vertex_ids)
      cell->centroid = cell->centroid + *umesh->vertices[vid];

    cell->centroid = cell->centroid/(cell->vertex_ids.size());
  }

  chi_log.Log(LOG_0) << "Computed centroids";
  MPI_Barrier(MPI_COMM_WORLD);


  //======================================== Apply partitioning scheme
  auto grid = chi_mesh::MeshContinuum::New();

  if (options.partition_type == PartitionType::KBA_STYLE_XY or
      options.partition_type == PartitionType::KBA_STYLE_XYZ)
    KBA(umesh, grid);
  else
    PARMETIS(umesh,grid);

  chi_log.Log(LOG_0) << "Cells loaded.";
  MPI_Barrier(MPI_COMM_WORLD);

  AddContinuumToRegion(grid, *mesh_handler->region_stack.back());


  chi_log.Log(LOG_0)
    << "VolumeMesherPredefined3D: Number of nodes in region = "
    << grid->vertices.size()
    << std::endl;
  grid->vertices.shrink_to_fit();

  chi_log.Log(LOG_ALLVERBOSE_1)
    << "### LOCATION[" << chi_mpi.location_id
    << "] amount of local cells="
    << grid->local_cell_glob_indices.size();

  int total_local_cells = grid->local_cells.size();
  int total_global_cells = 0;

  MPI_Allreduce(&total_local_cells,
                &total_global_cells,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD);

  chi_log.Log(LOG_0)
    << "VolumeMesherPredefined3D: Cells created = "
    << total_global_cells
    << std::endl;



}