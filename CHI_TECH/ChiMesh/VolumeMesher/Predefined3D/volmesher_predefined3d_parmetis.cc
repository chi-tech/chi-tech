#include "volmesher_predefined3d.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell_polyhedron.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

#include "petsc.h"

//###################################################################
/**Determines if a chi_mesh::UnpartitionedMesh::LightWeightCell is a
 * neighbor to the current partition for ParMETIS-style partitioning.
 * This method loops over the faces of the lightweight cell and
 * determines the partition-id of each the neighbors. If the neighbor
 * has a partition id equal to that of the current process then
 * it means this reference cell is a neighbor.*/
bool chi_mesh::VolumeMesherPredefined3D::
  IsRawCellNeighborToPartitionParmetis(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell,
    const std::vector<int>& cell_pids)
{
  auto handler = chi_mesh::GetCurrentHandler();

  auto umesh = handler->unpartitionedmesh_stack.back();

  bool is_neighbor = false;
  for (const auto& face : lwcell.faces)
  {
    if (face.neighbor < 0) continue;
    auto adj_cell = umesh->raw_cells[face.neighbor];
    int partition_id = cell_pids[face.neighbor];
    if (partition_id == chi_mpi.location_id)
    {
      is_neighbor = true;
      break;
    }
  }

  return is_neighbor;
}

//###################################################################
/** Applies KBA-style partitioning to the mesh.*/
void chi_mesh::VolumeMesherPredefined3D::
  PARMETIS(chi_mesh::UnpartitionedMesh* umesh,
           chi_mesh::MeshContinuum* grid)
{
  chi_log.Log(LOG_0) << "Partitioning mesh.";
  std::vector<int> cell_pids(umesh->raw_cells.size(),0);
  if (chi_mpi.location_id == 0)
  {
    //======================================== Build indices
    std::vector<int> i_indices(umesh->raw_cells.size()+1,0);
    std::vector<int> j_indices;
    int i=-1;
    int icount = 0;
    for (auto cell : umesh->raw_cells)
    {
      ++i;
      i_indices[i] = icount;

      for (auto& face : cell->faces)
        if (face.neighbor >= 0)
        {
          j_indices.push_back(face.neighbor);
          ++icount;
        }
    }
    i_indices[i+1] = icount;
    chi_log.Log(LOG_0VERBOSE_1) << "Done building indices.";

    //======================================== Copy to raw arrays
    int* i_indices_raw; // = new int[i_indices.size()];
    int* j_indices_raw; // = new int[j_indices.size()];
    PetscMalloc(i_indices.size()*sizeof(int),&i_indices_raw);
    PetscMalloc(j_indices.size()*sizeof(int),&j_indices_raw);

    for (int j=0; j<i_indices.size(); ++j)
      i_indices_raw[j] = i_indices[j];

    for (int j=0; j<j_indices.size(); ++j)
      j_indices_raw[j] = j_indices[j];

    chi_log.Log(LOG_0VERBOSE_1) << "Done copying to raw indices.";

    //========================================= Create adjacency matrix
    Mat Adj; //Adjacency matrix
    MatCreateMPIAdj(PETSC_COMM_SELF,
                    (int)umesh->raw_cells.size(),
                    (int)umesh->raw_cells.size(),
                    i_indices_raw,j_indices_raw,NULL,&Adj);

    chi_log.Log(LOG_0VERBOSE_1) << "Done creating adjacency matrix.";

    //========================================= Create partitioning
    MatPartitioning part;
    IS is,isg;
    MatPartitioningCreate(MPI_COMM_SELF,&part);
    MatPartitioningSetAdjacency(part,Adj);
    MatPartitioningSetType(part,"parmetis");
    MatPartitioningSetNParts(part,chi_mpi.process_count);
    MatPartitioningApply(part,&is);
    MatPartitioningDestroy(&part);
    MatDestroy(&Adj);
    ISPartitioningToNumbering(is,&isg);
    chi_log.Log(LOG_0VERBOSE_1) << "Done building paritioned index set.";

    //========================================= Get cell global indices
    const int* cell_pids_raw;
    ISGetIndices(is,&cell_pids_raw);
    i=0;
    for (auto cell : umesh->raw_cells)
    {
      cell_pids[i] = cell_pids_raw[i];
      ++i;
    }
    ISRestoreIndices(is,&cell_pids_raw);

    chi_log.Log(LOG_0VERBOSE_1) << "Done retrieving cell global indices.";
  }

  //======================================== Broadcast partitioning to all
  //                                         locations
  MPI_Bcast(cell_pids.data(),        //buffer [IN/OUT]
            umesh->raw_cells.size(), //count
            MPI_INT,                 //data type
            0,                       //root
            MPI_COMM_WORLD);         //communicator
  chi_log.Log(LOG_0) << "Done partitioning mesh.";

  //======================================== Load up the vertices
  for (auto vert : umesh->vertices)
    grid->vertices.push_back(new chi_mesh::Vertex(*vert));

  MPI_Barrier(MPI_COMM_WORLD);

  //======================================== Load up the cells
  int global_id=-1;
  for (auto raw_cell : umesh->raw_cells)
  {
    ++global_id;
    auto temp_cell = new chi_mesh::Cell(chi_mesh::CellType::GHOST);
    temp_cell->centroid = raw_cell->centroid;
    temp_cell->global_id = global_id;
    temp_cell->partition_id = cell_pids[global_id];
    temp_cell->material_id = raw_cell->material_id;

    if (temp_cell->partition_id != chi_mpi.location_id)
    {
      if (IsRawCellNeighborToPartitionParmetis(*raw_cell,cell_pids))
        grid->cells.push_back(temp_cell);
      else
        delete temp_cell;
    }
    else
    {
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
}