#include "volmesher_predefined3d.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/SurfaceMesher/surfacemesher.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"
#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell_polyhedron.h"
#include "ChiMesh/Region//chi_region.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog chi_log;
extern ChiMPI chi_mpi;

#include <ChiTimer/chi_timer.h>
extern ChiTimer chi_program_timer;

#include <ChiConsole/chi_console.h>
extern ChiConsole  chi_console;

//###################################################################
/**Gets the partition ID from a centroid.*/
int chi_mesh::VolumeMesherPredefined3D::
  GetPartitionIDFromCentroid(const chi_mesh::Vertex& centroid)
{
  auto handler = chi_mesh::GetCurrentHandler();

  int Px = handler->surface_mesher->partitioning_x;
  int Py = handler->surface_mesher->partitioning_y;


  chi_mesh::Cell temp_cell(chi_mesh::CellType::GHOST);
  temp_cell.centroid = centroid;

  auto xyz = GetCellXYZPartitionID(&temp_cell);

  int nxi = std::get<0>(xyz);
  int nyi = std::get<1>(xyz);
  int nzi = std::get<2>(xyz);

  return nzi*Px*Py + nyi*Px + nxi;
}

//###################################################################
/**Determines if a chi_mesh::UnpartitionedMesh::LightWeightCell is a
 * neighbor to the current partition.
 * This method loops over the faces of the lightweight cell and
 * determines the partition-id of each the neighbors. If the neighbor
 * has a partition id equal to that of the current process then
 * it means this reference cell is a neighbor.*/
bool chi_mesh::VolumeMesherPredefined3D::
  IsRawCellNeighborToPartition(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell)
{
  auto handler = chi_mesh::GetCurrentHandler();

  auto umesh = handler->unpartitionedmesh_stack.back();

  bool is_neighbor = false;
  for (const auto& face : lwcell.faces)
  {
    if (face.neighbor < 0) continue;
    auto adj_cell = umesh->raw_cells[face.neighbor];
    int partition_id = GetPartitionIDFromCentroid(adj_cell->centroid);
    if (partition_id == chi_mpi.location_id)
    {
      is_neighbor = true;
      break;
    }
  }

  return is_neighbor;
}

//###################################################################
/**Executes the predefined3D mesher.*/
void chi_mesh::VolumeMesherPredefined3D::Execute()
{
  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " VolumeMesherPredefined3D executed. Memory in use = "
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

      bool stop_searching = false;
      for (auto& adj_cell_id_set : cells_to_search)
      {
        for (auto& adj_cell_id : adj_cell_id_set)
        {
          if (adj_cell_id == c) continue;
          auto adj_cell = umesh->raw_cells[adj_cell_id];

          //Assume it matches now disprove
          bool adj_cell_matches = true;
          for (auto cfvid : face.vertex_ids)
          {
            bool vertex_found=false;
            for (auto acvid : adj_cell->vertex_ids)
              if (cfvid == acvid) {vertex_found = true; break;}

            if (not vertex_found) { adj_cell_matches = false; break;}
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


  //======================================== Load up the vertices
  auto grid = new chi_mesh::MeshContinuum;
  for (auto vert : umesh->vertices)
    grid->vertices.push_back(new chi_mesh::Vertex(*vert));

  chi_log.Log(LOG_0) << "Vertices loaded.";
  MPI_Barrier(MPI_COMM_WORLD);

  int loc_id = chi_mpi.location_id;

  //======================================== Load up the cells
  int global_id=-1;
  for (auto raw_cell : umesh->raw_cells)
  {
    ++global_id;
    auto temp_cell = new chi_mesh::Cell(chi_mesh::CellType::GHOST);
    temp_cell->centroid = raw_cell->centroid;
    temp_cell->global_id = global_id;
    temp_cell->partition_id = GetPartitionIDFromCentroid(temp_cell->centroid);
    temp_cell->material_id = raw_cell->material_id;

//    printf("[%d] Bla\n",loc_id);

    if (temp_cell->partition_id != chi_mpi.location_id)
    {
//      printf("[%d] Bla1\n",loc_id);
      if (IsRawCellNeighborToPartition(*raw_cell))
        grid->cells.push_back(temp_cell);
      else
        delete temp_cell;

//      printf("[%d] Bla2\n",loc_id);
    }
    else
    {
//      printf("[%d] Bla3\n",loc_id);
      auto polyh_cell = new chi_mesh::CellPolyhedron;
      polyh_cell->centroid = temp_cell->centroid;
      polyh_cell->global_id = temp_cell->global_id;
      polyh_cell->partition_id = temp_cell->partition_id;
      polyh_cell->material_id = temp_cell->material_id;

      polyh_cell->vertex_ids = raw_cell->vertex_ids;

      for (auto& raw_face : raw_cell->faces)
      {
        chi_mesh::CellFace newFace;

        newFace.neighbor = raw_face.neighbor;

        newFace.vertex_ids = raw_face.vertex_ids;
        auto vfc = chi_mesh::Vertex(0.0, 0.0, 0.0);
        for (auto fvid : newFace.vertex_ids)
          vfc = vfc + *grid->vertices[fvid];
        newFace.centroid = vfc / newFace.vertex_ids.size();

        //Compute normal
//        auto va = *grid->vertices[newFace.vertex_ids[0]] - vfc;
//        auto vb = *grid->vertices[newFace.vertex_ids[1]] - vfc;
//
//        auto vn = va.Cross(vb);
//        newFace.normal = vn/vn.Norm();

        newFace.normal = chi_mesh::Normal(0.0,0.0,0.0);
        int last_vert_ind = newFace.vertex_ids.size()-1;
        for (int fv=0; fv<newFace.vertex_ids.size(); ++fv)
        {
          int fvid_m = newFace.vertex_ids[fv];
          int fvid_p = (fv == last_vert_ind)? newFace.vertex_ids[0] :
                                              newFace.vertex_ids[fv+1];
          auto leg_m = *grid->vertices[fvid_m] - newFace.centroid;
          auto leg_p = *grid->vertices[fvid_p] - newFace.centroid;

          auto vn = leg_m.Cross(leg_p);

          newFace.normal = newFace.normal + vn.Normalized();
        }
        newFace.normal = (newFace.normal/newFace.vertex_ids.size()).Normalized();

        polyh_cell->faces.push_back(newFace);
      }

      grid->cells.push_back(polyh_cell);
    }//else
  }

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