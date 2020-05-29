#include "fv.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

//###################################################################
/**Builds finite volume based sparsity pattern.*/
void SpatialDiscretization_FV::BuildSparsityPattern(
  chi_mesh::MeshContinuum *grid,
  std::vector<int> &nodal_nnz_in_diag,
  std::vector<int> &nodal_nnz_off_diag,
  chi_math::UnknownManager* unknown_manager)
{
  unsigned int N = 1; //Number of components

  if (unknown_manager != nullptr)
    N = unknown_manager->GetTotalUnknownSize();

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  const size_t num_local_cells = grid->local_cells.size();
  block_size_per_unknown = num_local_cells;

  nodal_nnz_in_diag.resize(num_local_cells*N,0.0);
  nodal_nnz_off_diag.resize(num_local_cells*N,0.0);

  for (int block=0; block<N; ++block)
  {
    for (auto& cell : grid->local_cells)
    {
      int i=cell.local_id + block*num_local_cells;

      nodal_nnz_in_diag[i]   += 1;

      for (auto& face : cell.faces)
      {
        if (face.neighbor < 0) continue;

        if (face.IsNeighborLocal(grid))
          nodal_nnz_in_diag[i] += 1;
        else
          nodal_nnz_off_diag[i] += 1;
      }
    }//for cell
  }//for block

}