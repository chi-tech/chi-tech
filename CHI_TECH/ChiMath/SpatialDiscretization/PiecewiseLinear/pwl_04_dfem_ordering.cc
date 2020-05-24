#include "pwl.h"

#include "ChiTimer/chi_timer.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Reorders the nodes for parallel computation in a Continuous
 * Finite Element calculation.*/
std::pair<int,int> SpatialDiscretization_PWL::
  OrderNodesDFEM(chi_mesh::MeshContinuum *grid)
{
  ChiTimer t_stage[6];

  t_stage[0].Reset();
  //================================================== Get local DOF count
  size_t num_loc_cells = grid->local_cell_glob_indices.size();

  cell_dfem_block_address.resize(num_loc_cells,0);

  int local_dof_count=0;
  for (int lc=0; lc<num_loc_cells; lc++)
  {
    auto cell_fe_view = cell_fe_views[lc];
    cell_dfem_block_address[lc] = local_dof_count;
    local_dof_count += cell_fe_view->dofs;
  }

  //================================================== Get global DOF count
  int global_dof_count=0;
  MPI_Allreduce(&local_dof_count,    //Send buffer
                &global_dof_count,   //Recv buffer
                1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  //================================================== Ring communicate DOF start
  dfem_local_block_address = 0;
  if (chi_mpi.location_id != 0)
  {
    MPI_Recv(&dfem_local_block_address,
             1,MPI_INT,              //Count and type
             chi_mpi.location_id-1,  //Source
             111,                    //Tag
             MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }

  if (chi_mpi.location_id != (chi_mpi.process_count-1))
  {
    int next_loc_start = dfem_local_block_address + local_dof_count;
    MPI_Send(&next_loc_start,
             1,MPI_INT,
             chi_mpi.location_id+1,
             111,
             MPI_COMM_WORLD);
  }

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "Local dof count, start, total "
    << local_dof_count << " "
    << dfem_local_block_address << " "
    << global_dof_count;

  return {local_dof_count, global_dof_count};
}

//###################################################################
/**Provides a mapping of cell's DOF from a DFEM perspective.*/
int SpatialDiscretization_PWL::MapDFEMDOF(chi_mesh::Cell *cell,
                                          int dof,
                                          int component,
                                          int component_block_offset)
{
  if (cell->partition_id == chi_mpi.location_id)
  {
    int address = cell_dfem_block_address[cell->local_id] +
                  dfem_local_block_address +
                  dof;
    return address*component_block_offset + component;
  }
  else
  {
    int index = 0;
    bool found = false;
    for (auto neighbor_info : neighbor_cell_block_address)
    {
      if (neighbor_info.first == cell->global_id) {
        found = true; break;
      }
      ++index;
    }

    if (!found)
    {
      chi_log.Log(LOG_ALLERROR)
        << "SpatialDiscretization_PWL::MapDFEMDOF. Mapping failed for cell "
        << "with global index " << cell->global_id << " and partition-ID "
        << cell->partition_id;
      exit(EXIT_FAILURE);
    }

    int address = neighbor_cell_block_address[index].second + dof;
    return address*component_block_offset + component;
  }
}

//###################################################################
/**Provides a local mapping of cell's DOF from a DFEM perspective.*/
int SpatialDiscretization_PWL::MapDFEMDOFLocal(chi_mesh::Cell *cell,
                                               int dof,
                                               int component,
                                               int component_block_offset)
{
  if (cell->partition_id == chi_mpi.location_id)
  {
    int address = cell_dfem_block_address[cell->local_id] + dof;
    return address*component_block_offset + component;
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "SpatialDiscretization_PWL::MapDFEMDOF. Mapping failed for cell "
      << "with global index " << cell->global_id << " and partition-ID "
      << cell->partition_id;
    exit(EXIT_FAILURE);
  }
}