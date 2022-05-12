#include "fv.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;


//###################################################################
/**Develops node ordering per location.*/
void SpatialDiscretization_FV::
  OrderNodes()
{
  chi_log.Log(LOG_ALLVERBOSE_1) << "FV discretization - Reordering nodes.";

  fv_local_block_address = 0;

  if (chi::mpi.location_id != 0)
  {
    MPI_Recv(&fv_local_block_address,
             1,MPI_INT,              //Count and type
             chi::mpi.location_id-1,  //Source
             111,                    //Tag
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }

  if (chi::mpi.location_id != (chi::mpi.process_count-1))
  {
    int next_loc_start = fv_local_block_address + ref_grid->local_cells.size();
    MPI_Send(&next_loc_start,
             1,MPI_INT,
             chi::mpi.location_id+1,
             111,
             MPI_COMM_WORLD);
  }

  //======================================== Collect block addresses
  locJ_block_address.clear();
  locJ_block_address.resize(chi::mpi.process_count, 0);
  MPI_Allgather(&fv_local_block_address,      //send buf
                1,                            //send count
                MPI_INT,                      //send type
                locJ_block_address.data(), //recv buf
                1,                            //recv count
                MPI_INT,                      //recv type
                MPI_COMM_WORLD);              //communicator

  //======================================== Collect block sizes
  locJ_block_size.clear();
  locJ_block_size.resize(chi::mpi.process_count, 0);
  int num_local_cells =  ref_grid->local_cells.size();
  MPI_Allgather(&num_local_cells,             //send buf
                1,                            //send count
                MPI_INT,                      //send type
                locJ_block_size.data(),    //recv buf
                1,                            //recv count
                MPI_INT,                      //recv type
                MPI_COMM_WORLD);              //communicator

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "Local dof count, start "
    << ref_grid->local_cells.size() << " "
    << fv_local_block_address;
}

