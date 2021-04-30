#include "volmesher_predefunpart.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

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
bool chi_mesh::VolumeMesherPredefinedUnpartitioned::
  IsRawCellNeighborToPartitionParmetis(
    const chi_mesh::UnpartitionedMesh::LightWeightCell& lwcell,
    const std::vector<int64_t>& cell_pids)
{
  bool is_neighbor = false;
  for (const auto& face : lwcell.faces)
  {
    if (not face.neighbor) continue;
    int64_t partition_id = cell_pids[face.neighbor];
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
void chi_mesh::VolumeMesherPredefinedUnpartitioned::
  PARMETIS(chi_mesh::UnpartitionedMesh* umesh,
           chi_mesh::MeshContinuumPtr& grid)
{
  chi_log.Log(LOG_0) << "Partitioning mesh.";
  std::vector<int64_t> cell_pids(umesh->raw_cells.size(),0);
  if (chi_mpi.location_id == 0)
  {
    if (umesh->raw_cells.size() > 1)
    {
      //======================================== Build indices
      std::vector<int64_t> i_indices(umesh->raw_cells.size()+1,0);
      std::vector<int64_t> j_indices;
      int64_t i=-1;
      int64_t icount = 0;
      for (auto cell : umesh->raw_cells)
      {
        ++i;
        i_indices[i] = icount;

        for (auto& face : cell->faces)
          if (face.has_neighbor)
          {
            j_indices.push_back(int64_t(face.neighbor));
            ++icount;
          }
      }
      i_indices[i+1] = icount;
      chi_log.Log(LOG_0VERBOSE_1) << "Done building indices.";

      //======================================== Copy to raw arrays
      int64_t* i_indices_raw;
      int64_t* j_indices_raw;
      PetscMalloc(i_indices.size()*sizeof(int64_t),&i_indices_raw);
      PetscMalloc(j_indices.size()*sizeof(int64_t),&j_indices_raw);

      for (int64_t j=0; j<i_indices.size(); ++j)
        i_indices_raw[j] = i_indices[j];

      for (int64_t j=0; j<j_indices.size(); ++j)
        j_indices_raw[j] = j_indices[j];

      chi_log.Log(LOG_0VERBOSE_1) << "Done copying to raw indices.";

      //========================================= Create adjacency matrix
      Mat Adj; //Adjacency matrix
      MatCreateMPIAdj(PETSC_COMM_SELF,
                      (int64_t)umesh->raw_cells.size(),
                      (int64_t)umesh->raw_cells.size(),
                      i_indices_raw, j_indices_raw, nullptr, &Adj);

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
      const int64_t* cell_pids_raw;
      ISGetIndices(is,&cell_pids_raw);
      i=0;
      for (auto cell : umesh->raw_cells)
      {
        cell_pids[i] = cell_pids_raw[i];
        ++i;
      }
      ISRestoreIndices(is,&cell_pids_raw);

      chi_log.Log(LOG_0VERBOSE_1) << "Done retrieving cell global indices.";
    }//if more than 1 cell

  }

  //======================================== Broadcast partitioning to all
  //                                         locations
  MPI_Bcast(cell_pids.data(),        //buffer [IN/OUT]
            umesh->raw_cells.size(), //count
            MPI_LONG_LONG_INT,       //data type
            0,                       //root
            MPI_COMM_WORLD);         //communicator
  chi_log.Log(LOG_0) << "Done partitioning mesh.";

  //======================================== Load up the vertices
  for (auto& vert : umesh->vertices)
    grid->vertices.push_back(vert);

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
      if (raw_cell->type == chi_mesh::CellType::SLAB)
        AddSlabToGrid(*raw_cell,*temp_cell,*grid);
      else if (raw_cell->type == chi_mesh::CellType::POLYGON)
        AddPolygonToGrid(*raw_cell,*temp_cell,*grid);
      else if (raw_cell->type == chi_mesh::CellType::POLYHEDRON)
        AddPolyhedronToGrid(*raw_cell,*temp_cell,*grid);
      else
      {
        chi_log.Log(LOG_ALLERROR)
          << "Unsupported cell encountered in "
             "chi_mesh::VolumeMesherPredefinedUnpartitioned";
        exit(EXIT_FAILURE);
      }
    }//else
  }
}