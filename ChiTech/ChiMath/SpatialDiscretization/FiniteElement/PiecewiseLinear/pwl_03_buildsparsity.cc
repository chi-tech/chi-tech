#include "pwl.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Builds the sparsity pattern for a Discontinuous Finite Element Method.*/
void SpatialDiscretization_PWLD::
BuildSparsityPattern(std::vector<int64_t> &nodal_nnz_in_diag,
                     std::vector<int64_t> &nodal_nnz_off_diag,
                     chi_math::UnknownManager& unknown_manager)
{
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL CONNECTIVITY
  size_t local_dof_count = local_base_block_size;

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag .resize(local_dof_count,0);
  nodal_nnz_off_diag.resize(local_dof_count,0);

  int lc=0;
  for (auto& cell : ref_grid->local_cells)
  {
    auto cell_fe_view = GetCellMappingFE(lc);
    size_t num_nodes = cell_fe_view->num_nodes;

    //==================================== Self connection
    for (int i=0; i<num_nodes; ++i)
    {
      int ir = cell_local_block_address[lc] + i;
      nodal_nnz_in_diag[ir] += num_nodes;
    }

    //==================================== Local adjacent cell connections
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor and face.IsNeighborLocal(*ref_grid))
      {
        auto& adj_cell = ref_grid->cells[face.neighbor_id];
        auto adj_cell_fe_view = GetCellMappingFE(adj_cell.local_id);

        for (int i=0; i<num_nodes; ++i)
        {
          int ir = cell_local_block_address[lc] + i;
          nodal_nnz_in_diag[ir] += adj_cell_fe_view->num_nodes;
        }
      }
    }//for face
    ++lc;
  }//for local cell

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORING CONNECTIVITY
  lc=0;
  for (auto& cell : ref_grid->local_cells)
  {
    auto cell_fe_view = GetCellMappingFE(lc);

    //==================================== Local adjacent cell connections
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor and (not face.IsNeighborLocal(*ref_grid)))
      {
        auto adj_cell_fe_view = GetNeighborCellMappingFE(face.neighbor_id);

        for (int i=0; i<cell_fe_view->num_nodes; ++i)
        {
          int ir = cell_local_block_address[lc] + i;
          nodal_nnz_off_diag[ir] += adj_cell_fe_view->num_nodes;
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
  chi_log.Log(LOG_0) << "Done building DFEM sparsity pattern";

}





