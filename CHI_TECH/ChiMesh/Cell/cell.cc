#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include "cell.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog chi_log;
extern ChiMPI chi_mpi;

//###################################################################
/**Determines the neighbor's partition and whether its local or not.*/
bool chi_mesh::CellFace::
  IsNeighborLocal(chi_mesh::MeshContinuum *grid)
{
  if (neighbor < 0) return false;
  if (chi_mpi.process_count == 1) return true;

  if (not neighbor_parallel_info_initialized)
    InitializeNeighborParallelInfo(grid);


  return (neighbor_partition_id == chi_mpi.location_id);
}

//###################################################################
/**Determines the neighbor's partition.*/
int chi_mesh::CellFace::
  GetNeighborPartitionID(chi_mesh::MeshContinuum *grid)
{
  if (neighbor < 0) return -1;
  if (chi_mpi.process_count == 1) return 0;

  if (not neighbor_parallel_info_initialized)
    InitializeNeighborParallelInfo(grid);

  return neighbor_partition_id;
}

//###################################################################
/**Determines the neighbor's local id.*/
int chi_mesh::CellFace::
  GetNeighborLocalID(chi_mesh::MeshContinuum *grid)
{
  if (neighbor < 0) return -1;
  if (chi_mpi.process_count == 1) return neighbor;

  if (not neighbor_parallel_info_initialized)
    InitializeNeighborParallelInfo(grid);

  return neighbor_local_id;
}

//###################################################################
/**Determines the neighbor's associated face.*/
int chi_mesh::CellFace::
  GetNeighborAssociatedFace(chi_mesh::MeshContinuum *grid)
{
  auto& cur_face = *this;
  //======================================== Check index validity
  if ((cur_face.neighbor<0) || (not cur_face.IsNeighborLocal(grid)))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid cell index encountered in call to "
      << "CellFace::GetNeighborAssociatedFace. Index points "
      << "to either a boundary"
      << "or a non-local cell.";
    exit(EXIT_FAILURE);
  }

  auto adj_cell = &grid->local_cells[cur_face.GetNeighborLocalID(grid)];

  int associated_face = -1;
  if (cur_face.neighbor_ass_face < 0)
  {
    //======================================== Loop over adj cell faces
    int af=-1;
    for (auto& adj_face : adj_cell->faces)
    {
      ++af;
      //Assume face matches
      bool face_matches = true; //Now disprove it
      //================================= Loop over adj cell face verts
      for (auto afvi : adj_face.vertex_ids)
      {
        //========================== Try and find them in the reference face
        bool found = false;
        for (auto cfvi : cur_face.vertex_ids)
        {
          if (cfvi == afvi)
          {
            found = true;
            break;
          }
        }//for cfv

        if (!found) {face_matches = false; break;}
      }//for afv

      if (face_matches) {associated_face = af; break;}
    }
  }
  else
    associated_face = cur_face.neighbor_ass_face;

  //======================================== Check associated face validity
  if (associated_face<0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Could not find associated face in call to "
      << "CellFace::GetNeighborAssociatedFace. Reference face with centroid at \n"
      << cur_face.centroid.PrintS();
    for (int af=0; af < adj_cell->faces.size(); af++)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Adjacent cell face " << af << " centroid "
        << adj_cell->faces[af].centroid.PrintS();
    }
    exit(EXIT_FAILURE);
  }

  cur_face.neighbor_ass_face = associated_face;
  return associated_face;
}

//###################################################################
/**Initializes a neighbor's parallel information.*/
void chi_mesh::CellFace::
  InitializeNeighborParallelInfo(chi_mesh::MeshContinuum *grid)
{
  auto adj_cell = grid->cells[neighbor];
  neighbor_partition_id = adj_cell->partition_id;
  neighbor_parallel_info_initialized = true;

  if (neighbor_partition_id == chi_mpi.location_id)
    neighbor_local_id = adj_cell->local_id;
}