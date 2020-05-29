#include "fv.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Develops node ordering per location.*/
void SpatialDiscretization_FV::
  ReOrderNodes(chi_mesh::MeshContinuum* grid)
{
  chi_log.Log(LOG_ALLVERBOSE_1) << "FV discretization - Reordering nodes.";

  fv_local_block_address = 0;

  if (chi_mpi.location_id != 0)
  {
    MPI_Recv(&fv_local_block_address,
             1,MPI_INT,              //Count and type
             chi_mpi.location_id-1,  //Source
             111,                    //Tag
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }

  if (chi_mpi.location_id != (chi_mpi.process_count-1))
  {
    int next_loc_start = fv_local_block_address + grid->local_cells.size();
    MPI_Send(&next_loc_start,
             1,MPI_INT,
             chi_mpi.location_id+1,
             111,
             MPI_COMM_WORLD);
  }

  //======================================== Collect block addresses
  locJ_fv_block_address.clear();
  locJ_fv_block_address.resize(chi_mpi.process_count,0);
  MPI_Allgather(&fv_local_block_address,      //send buf
                1,                            //send count
                MPI_INT,                      //send type
                locJ_fv_block_address.data(), //recv buf
                1,                            //recv count
                MPI_INT,                      //recv type
                MPI_COMM_WORLD);              //communicator

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "Local dof count, start "
    << grid->local_cells.size() << " "
    << fv_local_block_address;
}

//###################################################################
/**Maps a finite volume degree of freedom.*/
int SpatialDiscretization_FV::
MapDOF(chi_mesh::Cell* cell,
       chi_math::UnknownManager* unknown_manager,
       unsigned int unknown_id,
       unsigned int component)
{
  if (cell == nullptr)
  {
    chi_log.Log(LOG_ALLERROR)
      << "SpatialDiscretization_FV::MapDOF reference cell is nullptr.";
    exit(EXIT_FAILURE);
  }
  if (component < 0) return -1;

  size_t num_unknowns = unknown_manager->GetTotalUnknownSize();

  if (component >= num_unknowns) return -1;

  size_t block_id = unknown_manager->MapUnknown(unknown_id,component);

  int address=-1;
  if (cell->partition_id == chi_mpi.location_id)
    address = fv_local_block_address + cell->local_id;
  else
  {
//    int ghost_local_id = ref_grid->cells.GetGhostLocalID(cell->global_id);
//    auto ghost_cell = neighbor_cells[ghost_local_id];
//
//    address = locJ_fv_block_address[ghost_cell->partition_id] +
//              ghost_cell->local_id;
    int ghost_local_id = 0;
    for (auto ghost : neighbor_cells)
      if (ghost->global_id == cell->global_id)
        ghost_local_id = ghost->local_id;

    address = locJ_fv_block_address[cell->partition_id] +
              ghost_local_id;
  }

  return num_unknowns*address + block_id;
}

//###################################################################
/**Maps a finite volume degree of freedom to a local address.*/
int SpatialDiscretization_FV::
MapDOFLocal(chi_mesh::Cell* cell,
            chi_math::UnknownManager* unknown_manager,
            unsigned int unknown_id,
            unsigned int component)
{
  if (component < 0) return -1;

  size_t num_unknowns = unknown_manager->GetTotalUnknownSize();

  if (component >= num_unknowns) return -1;

  size_t block_id = unknown_manager->MapUnknown(unknown_id,component);

  int address = -1;
  if (cell->partition_id == chi_mpi.location_id)
    address = cell->local_id;
  else
    address = ref_grid->local_cells.size() +
              ref_grid->cells.GetGhostLocalID(cell->global_id);

  if (address >=0)
    return num_unknowns*address + block_id;
  else
    return -1;
}