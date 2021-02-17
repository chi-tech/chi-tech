#include "pwl.h"

#include "ChiTimer/chi_timer.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Reorders the nodes for parallel computation in a Continuous
 * Finite Element calculation.*/
void SpatialDiscretization_PWL::
  OrderNodes(chi_mesh::MeshContinuumPtr grid)
{
  ChiTimer t_stage[6];

  t_stage[0].Reset();
  //================================================== Check cell views avail
  size_t num_loc_cells = grid->local_cell_glob_indices.size();
  if (num_loc_cells != cell_fe_views.size())
  {
    chi_log.Log(LOG_ALLERROR)
      << "SpatialDiscretization_PWL::OrderNodesDFEM. Number of cell_fe_views ("
      << cell_fe_views.size() << ") does not correspond to the number of cells ("
      << num_loc_cells << "). Was a call to AddViewOfLocalContinuum made?";
    exit(EXIT_FAILURE);
  }


  //================================================== Get local DOF count
  cell_local_block_address.resize(num_loc_cells, 0);

  int local_dof_count=0;
  for (int lc=0; lc<num_loc_cells; lc++)
  {
    auto cell_fe_view = cell_fe_views[lc];
    cell_local_block_address[lc] = local_dof_count;
    local_dof_count += cell_fe_view->dofs;
  }

  //================================================== Get global DOF count
  int global_dof_count=0;
  MPI_Allreduce(&local_dof_count,    //Send buffer
                &global_dof_count,   //Recv buffer
                1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  //================================================== Ring communicate DOF start
  local_block_address = 0;
  if (chi_mpi.location_id != 0)
  {
    MPI_Recv(&local_block_address,
             1,MPI_INT,              //Count and type
             chi_mpi.location_id-1,  //Source
             111,                    //Tag
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }

  if (chi_mpi.location_id != (chi_mpi.process_count-1))
  {
    int next_loc_start = local_block_address + local_dof_count;
    MPI_Send(&next_loc_start,
             1,MPI_INT,
             chi_mpi.location_id+1,
             111,
             MPI_COMM_WORLD);
  }

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "Local dof count, start, total "
    << local_dof_count << " "
    << local_block_address << " "
    << global_dof_count;

  local_base_block_size = local_dof_count;
  globl_base_block_size = global_dof_count;

//  //======================================== Collect block addresses
//  locJ_block_address.clear();
//  locJ_block_address.resize(chi_mpi.process_count, 0);
//  MPI_Allgather(&dfem_local_block_address,    //send buf
//                1,                            //send count
//                MPI_INT,                      //send type
//                locJ_block_address.data(),    //recv buf
//                1,                            //recv count
//                MPI_INT,                      //recv type
//                MPI_COMM_WORLD);              //communicator

  //======================================== Collect block sizes
  locJ_block_size.clear();
  locJ_block_size.resize(chi_mpi.process_count, 0);
  MPI_Allgather(&local_base_block_size,       //send buf
                1,                            //send count
                MPI_INT,                      //send type
                locJ_block_size.data(),       //recv buf
                1,                            //recv count
                MPI_INT,                      //recv type
                MPI_COMM_WORLD);              //communicator
}

