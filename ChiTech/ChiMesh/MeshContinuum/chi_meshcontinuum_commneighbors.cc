#include "chi_meshcontinuum.h"

#include "ChiMesh/Cell/cell.h"

#include "ChiDataTypes/byte_array.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Communicates neighboring cells to this location for use by methods
 * such as the interior penalty method. The method populates the
 * supplied vector neighbor_cells. The complete cell is
 * not populated, the face neighbors are not transferred.*/
void chi_mesh::MeshContinuum::CommunicatePartitionNeighborCells(
  std::map<uint64_t, std::unique_ptr<chi_mesh::Cell>>& neighbor_cells)
{
  MPI_Barrier(MPI_COMM_WORLD);

  const auto& ghost_cell_ids = cells.GetGhostGlobalIDs();

  for (uint64_t global_id : ghost_cell_ids)
    neighbor_cells.insert(
      std::make_pair(global_id, &cells[global_id])
      );

}