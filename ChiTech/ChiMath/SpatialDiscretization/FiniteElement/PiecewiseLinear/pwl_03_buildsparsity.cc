#include "pwl.h"

#include <chi_log.h>
#include <chi_mpi.h>


//###################################################################
/**Builds the sparsity pattern for a Discontinuous Finite Element Method.*/
void chi_math::SpatialDiscretization_PWLD::
BuildSparsityPattern(std::vector<int64_t> &nodal_nnz_in_diag,
                     std::vector<int64_t> &nodal_nnz_off_diag,
                     const chi_math::UnknownManager& unknown_manager) const
{
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL CONNECTIVITY
  size_t local_dof_count = local_base_block_size;

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag .resize(local_dof_count,0);
  nodal_nnz_off_diag.resize(local_dof_count,0);

  int lc=0;
  for (const auto& cell : ref_grid->local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);
    size_t num_nodes = cell_mapping.NumNodes();

    //==================================== Self connection
    for (int i=0; i<num_nodes; ++i)
    {
      int64_t ir = cell_local_block_address[lc] + i;
      nodal_nnz_in_diag[ir] += static_cast<int64_t>(num_nodes);
    }

    //==================================== Local adjacent cell connections
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor and face.IsNeighborLocal(*ref_grid))
      {
        const auto& adj_cell = ref_grid->cells[face.neighbor_id];
        const auto& adj_cell_mapping = GetCellMapping(cell);

        for (int i=0; i<num_nodes; ++i)
        {
          int64_t ir = cell_local_block_address[lc] + i;
          nodal_nnz_in_diag[ir] +=
            static_cast<int64_t>(adj_cell_mapping.NumNodes());
        }
      }
    }//for face
    ++lc;
  }//for local cell

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORING CONNECTIVITY
  lc=0;
  for (auto& cell : ref_grid->local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);

    //==================================== Local adjacent cell connections
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor and (not face.IsNeighborLocal(*ref_grid)))
      {
        const auto& adj_cell = ref_grid->cells[face.neighbor_id];
        const auto& adj_cell_mapping = GetCellMapping(cell);

        for (int i=0; i < cell_mapping.NumNodes(); ++i)
        {
          int64_t ir = cell_local_block_address[lc] + i;
          nodal_nnz_off_diag[ir] +=
            static_cast<int64_t>(adj_cell_mapping.NumNodes());
        }
      }
    }
    ++lc;
  }//for local cell

  //============================================= Spacing according to unknown
  //                                              manager
  auto backup_nnz_in_diag  = nodal_nnz_in_diag;
  auto backup_nnz_off_diag = nodal_nnz_off_diag;

  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag.resize(local_base_block_size*N,0);
  nodal_nnz_off_diag.resize(local_base_block_size*N,0);

  if (unknown_manager.dof_storage_type == chi_math::UnknownStorageType::NODAL)
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
  else if (unknown_manager.dof_storage_type == chi_math::UnknownStorageType::BLOCK)
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

  MPI_Barrier(MPI_COMM_WORLD);
  chi::log.Log() << "Done building DFEM sparsity pattern";

}





