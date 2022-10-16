#include "fv.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

//###################################################################
/**Builds finite volume based sparsity pattern.*/
void chi_math::SpatialDiscretization_FV::BuildSparsityPattern(
  std::vector<int64_t> &nodal_nnz_in_diag,
  std::vector<int64_t> &nodal_nnz_off_diag,
  chi_math::UnknownManager& unknown_manager)
{
  unsigned int num_uk = unknown_manager.unknowns.size(); //Number of unknowns
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize(); //Total number of unknowns

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  const size_t num_local_cells = ref_grid->local_cells.size();

  nodal_nnz_in_diag.resize(num_local_cells*N,0.0);
  nodal_nnz_off_diag.resize(num_local_cells*N,0.0);

  for (int uk=0; uk<num_uk; ++uk)
  {
    const unsigned int num_comps = unknown_manager.unknowns[uk].num_components;
    for (int comp=0; comp<num_comps; ++comp)
    {
      for (auto& cell : ref_grid->local_cells)
      {
        const int64_t i = MapDOFLocal(cell,0,unknown_manager,uk,comp);

        nodal_nnz_in_diag[i]   += 1;

        for (auto& face : cell.faces)
        {
          if (not face.has_neighbor) continue;

          if (face.IsNeighborLocal(*ref_grid))
            nodal_nnz_in_diag[i] += 1;
          else
            nodal_nnz_off_diag[i] += 1;
        }
      }//for cell
    }//for components
  }//for unknown

}