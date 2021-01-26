#include "pwl.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Builds the sparsity pattern for a Discontinuous Finite Element Method.*/
void SpatialDiscretization_PWL::
  BuildDFEMSparsityPattern(chi_mesh::MeshContinuum *grid,
                           std::vector<int> &nodal_nnz_in_diag,
                           std::vector<int> &nodal_nnz_off_diag,
                           chi_math::UnknownManager* unknown_manager)
{
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL CONNECTIVITY
  int local_dof_count = local_base_block_size;

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag .resize(local_dof_count,0);
  nodal_nnz_off_diag.resize(local_dof_count,0);

  int lc=0;
  for (auto& cell : grid->local_cells)
  {
    auto cell_fe_view = cell_fe_views[lc];

    //==================================== Self connection
    for (int i=0; i<cell_fe_view->dofs; ++i)
    {
      int ir = cell_dfem_block_address[lc] + i;
      nodal_nnz_in_diag[ir] += cell_fe_view->dofs;
    }

    //==================================== Local adjacent cell connections
    for (auto& face : cell.faces)
    {
      if (face.IsNeighborLocal(grid))
      {
        int  adj_cell_local_id = face.GetNeighborLocalID(grid);
        auto adj_cell_fe_view = cell_fe_views[adj_cell_local_id];

        for (int i=0; i<cell_fe_view->dofs; ++i)
        {
          int ir = cell_dfem_block_address[lc] + i;
          nodal_nnz_in_diag[ir] += adj_cell_fe_view->dofs;
        }
      }
    }
    ++lc;
  }//for local cell



  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORING CONNECTIVITY
  lc=0;
  for (auto& cell : grid->local_cells)
  {
    auto cell_fe_view = cell_fe_views[lc];

    //==================================== Local adjacent cell connections
    for (auto& face : cell.faces)
    {
      if ((face.has_neighbor) and (not face.IsNeighborLocal(grid)) )
      {
        auto adj_cell_fe_view = MapNeighborCellFeView(face.neighbor_id);

        for (int i=0; i<cell_fe_view->dofs; ++i)
        {
          int ir = cell_dfem_block_address[lc] + i;
          nodal_nnz_off_diag[ir] += adj_cell_fe_view->dofs;
        }
      }
    }
    ++lc;
  }//for local cell



  std::set<int> local_neighboring_cell_indices;
  std::set<int> neighboring_partitions;

  //============================================= Collect local neighboring
  //                                              cell indices and neighboring
  //                                              partition indices
  // First establish which local nodes are
  // neighbors to other partitions. Build a set
  // of local cells and a set of destination
  // partitions.
  for (auto& cell : grid->local_cells)
  {
    for (auto& face : cell.faces)
    {
      if ((face.has_neighbor) and (not face.IsNeighborLocal(grid)) )
      {
        local_neighboring_cell_indices.insert(cell.local_id);
        neighboring_partitions.insert(face.GetNeighborPartitionID(grid));
      }
    }
  }

  //============================================= Subscribe neighbor-partitions
  //                                              the local cells they need
  // Now we establish a pair for each unique location.
  // pair.first is the partition-id of the location.
  // pair.second is the list of cells to be sent to that location.
  size_t num_neighbor_partitions = neighboring_partitions.size();

  typedef std::pair<int,std::vector<int>> ListOfCells;

  std::vector<ListOfCells> destination_subscriptions;
  destination_subscriptions.reserve(num_neighbor_partitions);
  for (int adj_part : neighboring_partitions)
  {
    ListOfCells new_list(-1,std::vector<int>());

    new_list.first = adj_part;
    for (int local_cell_index : local_neighboring_cell_indices)
    {
      auto& cell = grid->local_cells[local_cell_index];

      for (auto& face : cell.faces)
      {
        if ((face.has_neighbor) and (not face.IsNeighborLocal(grid)) )
        {
          if (face.GetNeighborPartitionID(grid) == adj_part)
            new_list.second.push_back(local_cell_index);
        }
      }//for faces
    }//for neighbor cells

    destination_subscriptions.push_back(new_list);
  }//for adj partition

  //============================================= Serialize
  // For each location we now serialize the cells
  // it needs.
  // The serialized values will be as follows
  // - cell_glob_index
  // - cell_block_address
  std::vector<ListOfCells> destination_serialized_data;

  destination_serialized_data.reserve(num_neighbor_partitions);

  int total_serial_size = 0;
  for (auto& cell_list : destination_subscriptions)
  {
    ListOfCells new_serial_data(-1,std::vector<int>());

    new_serial_data.first = cell_list.first;
    for (int local_cell_index : cell_list.second)
    {
      auto& cell = grid->local_cells[local_cell_index];

      std::vector<int>& border_cell_info = new_serial_data.second;

      int cell_global_block_address =
        cell_dfem_block_address[local_cell_index];

      border_cell_info.push_back(cell.global_id);         //cell_glob_index
      border_cell_info.push_back(cell_global_block_address);   //block address
    }
    total_serial_size += new_serial_data.second.size();
    destination_serialized_data.push_back(new_serial_data);
  }

  //============================================= Build send-counts and
  //                                              send-displacements arrays, and
  //                                              serialized vector
  std::vector<int> send_counts(chi_mpi.process_count,0);
  std::vector<int> send_displs(chi_mpi.process_count,0);
  std::vector<int> global_serialized_data;

  global_serialized_data.reserve(total_serial_size);

  int displacement = 0;
  for (auto& info : destination_serialized_data)
  {
    send_counts[info.first] = info.second.size();
    send_displs[info.first] = displacement;
    displacement += info.second.size();

    for (int val : info.second)
      global_serialized_data.push_back(val);
  }

  //============================================= Communicate counts
  std::vector<int> recv_counts(chi_mpi.process_count,0);

  MPI_Alltoall(send_counts.data(), 1, MPI_INT,
               recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

  //============================================= Build receive displacements
  std::vector<int> recv_displs(chi_mpi.process_count,0);
  int total_receive_size = 0;
  int c=0;
  for (auto val : recv_counts)
  {
    recv_displs[c] = total_receive_size;
    total_receive_size += val;
    ++c;
  }

  //============================================= Receive serialized data
  std::vector<int> global_receive_data(total_receive_size,0);
  MPI_Alltoallv(global_serialized_data.data(),
                send_counts.data(),
                send_displs.data(),
                MPI_INT,
                global_receive_data.data(),
                recv_counts.data(),
                recv_displs.data(),
                MPI_INT,
                MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Deserialize
  {
    int k=0;
    while (k<global_receive_data.size())
    {
      int cell_global_id     = global_receive_data[k]; k++;
      int cell_block_address = global_receive_data[k]; k++;

      neighbor_cell_block_address.emplace_back(cell_global_id, cell_block_address);
    }
  }

  //============================================= Spacing according to unknown
  //                                              manager
  if (unknown_manager != nullptr)
  {
    auto backup_nnz_in_diag  = nodal_nnz_in_diag;
    auto backup_nnz_off_diag = nodal_nnz_off_diag;

    unsigned int N = unknown_manager->GetTotalUnknownSize();

    nodal_nnz_in_diag.clear();
    nodal_nnz_off_diag.clear();

    nodal_nnz_in_diag.resize(local_base_block_size*N,0);
    nodal_nnz_off_diag.resize(local_base_block_size*N,0);

    if (unknown_manager->dof_storage_type == chi_math::DOFStorageType::NODAL)
    {
      int ir = -1;
      for (int i=0; i<local_base_block_size; ++i)
      {
        for (int j=0; j<N; ++j)
        {
          ++ir;
          nodal_nnz_in_diag[ir] = backup_nnz_in_diag[i];
          nodal_nnz_off_diag[ir] = backup_nnz_off_diag[i];
        }//for j
      }//for i
    }
    else if (unknown_manager->dof_storage_type == chi_math::DOFStorageType::BLOCK)
    {
      int ir = -1;
      for (int j=0; j<N; ++j)
      {
        for (int i=0; i<local_base_block_size; ++i)
        {
          ++ir;
          nodal_nnz_in_diag[ir] = backup_nnz_in_diag[i];
          nodal_nnz_off_diag[ir] = backup_nnz_off_diag[i];
        }//for i
      }//for j
    }
  }//if unknown manager supplied

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0) << "Done building DFEM sparsity pattern";

}