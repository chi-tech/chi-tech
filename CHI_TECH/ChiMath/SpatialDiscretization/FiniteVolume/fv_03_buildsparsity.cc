#include "fv.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

//###################################################################
/**Builds finite volume based sparsity pattern.*/
void SpatialDiscretization_FV::BuildSparsityPattern(
  chi_mesh::MeshContinuum* grid,
  std::vector<int> &nodal_nnz_in_diag,
  std::vector<int> &nodal_nnz_off_diag,
  chi_math::UnknownManager* unknown_manager)
{
  unsigned int num_uk = 1; //Number of unknowns
  unsigned int N = 1;      //Total number of unknowns

  if (unknown_manager != nullptr)
  {
    num_uk = unknown_manager->unknowns.size();
    N = unknown_manager->GetTotalUnknownSize();
  }

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  const size_t num_local_cells = grid->local_cells.size();

  nodal_nnz_in_diag.resize(num_local_cells*N,0.0);
  nodal_nnz_off_diag.resize(num_local_cells*N,0.0);

  for (int uk=0; uk<num_uk; ++uk)
  {
    int num_comps = unknown_manager->unknowns[uk].num_components;
    for (int comp=0; comp<num_comps; ++comp)
    {
      for (auto& cell : grid->local_cells)
      {
        int i = cell.local_id*N + comp;
        if (unknown_manager != nullptr)
          i=MapDOFLocal(&cell,unknown_manager,uk,comp);

        nodal_nnz_in_diag[i]   += 1;

        for (auto& face : cell.faces)
        {
          if (not face.has_neighbor) continue;

          if (face.IsNeighborLocal(*grid))
            nodal_nnz_in_diag[i] += 1;
          else
            nodal_nnz_off_diag[i] += 1;
        }
      }//for cell
    }//for components
  }//for unknown

}