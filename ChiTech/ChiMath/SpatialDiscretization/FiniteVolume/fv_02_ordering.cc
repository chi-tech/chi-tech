#include "fv.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
#include "chi_mpi.h"

//###################################################################
/**Develops node ordering per location.*/
void chi_math::SpatialDiscretization_FV::
  OrderNodes()
{
  chi::log.LogAllVerbose1() << "FV discretization - Reordering nodes.";

  //============================================= Communicate node counts
  const uint64_t local_num_nodes = ref_grid->local_cells.size();
  locJ_block_size.assign(chi::mpi.process_count,0);
  MPI_Allgather(&local_num_nodes,       // sendbuf
                1, MPI_UINT64_T,        // sendcount, sendtype
                locJ_block_size.data(), // recvbuf
                1, MPI_UINT64_T,        // recvcount, recvtype
                MPI_COMM_WORLD);        // comm

  //============================================= Build block addresses
  locJ_block_address.assign(chi::mpi.process_count,0);
  uint64_t global_num_nodes = 0;
  for (int j=0; j<chi::mpi.process_count; ++j)
  {
    locJ_block_address[j] = global_num_nodes;
    global_num_nodes += locJ_block_size[j];
  }

  local_block_address = locJ_block_address[chi::mpi.location_id];

  local_base_block_size = local_num_nodes;
  globl_base_block_size = global_num_nodes;


  chi::log.LogAllVerbose2()
    << "Local dof count, start "
    << ref_grid->local_cells.size() << " "
    << local_block_address;
}

