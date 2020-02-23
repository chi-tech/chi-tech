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

  if (not neighbor_partition_id_updated)
  {
    auto adj_cell = grid->cells[neighbor];
    neighbor_partition_id = adj_cell->partition_id;
    neighbor_partition_id_updated = true;
  }

  return (neighbor_partition_id == chi_mpi.location_id);
}

//###################################################################
/**Determines the neighbor's partition.*/
int chi_mesh::CellFace::
  GetNeighborPartitionID(chi_mesh::MeshContinuum *grid)
{
  if (neighbor < 0) return -1;
  if (chi_mpi.process_count == 1) return 0;

  if (not neighbor_partition_id_updated)
  {
    auto adj_cell = grid->cells[neighbor];
    neighbor_partition_id = adj_cell->partition_id;
    neighbor_partition_id_updated = true;
  }

  return neighbor_partition_id;
}

//###################################################################
/**Determines the neighbor's local id.*/
int chi_mesh::CellFace::
GetNeighborLocalID(chi_mesh::MeshContinuum *grid)
{
  if (neighbor < 0) return -1;
  if (chi_mpi.process_count == 1) return neighbor;

  if (not neighbor_local_id_updated)
  {
    if (IsNeighborLocal(grid))
    {
      auto adj_cell = grid->cells[neighbor];
      neighbor_local_id = adj_cell->cell_local_id;
    }

    neighbor_local_id_updated = true;
  }

  return neighbor_local_id;
}