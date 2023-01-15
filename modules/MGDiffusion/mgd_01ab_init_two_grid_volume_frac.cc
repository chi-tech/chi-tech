#include "mg_diffusion_solver.h"
#include "ChiTimer/chi_timer.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

//============================================= assemble matrix A
void mg_diffusion::Solver::Compute_TwoGrid_VolumeFractions()
{
  const auto& grid = *grid_ptr;
  const auto& sdm  = *sdm_ptr;

  const size_t ncells = grid.local_cells.size();
  VF.resize(ncells);

  int counter = 0;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto  qp_data      = cell_mapping.MakeVolumeQuadraturePointData();
    const size_t num_nodes   = cell_mapping.NumNodes();

    VF[counter].resize(num_nodes, 0.0);

    for (size_t i=0; i < num_nodes; ++i)
    {
      double vol_frac_shape_i = 0.0;
      for (size_t qp : qp_data.QuadraturePointIndices())
        vol_frac_shape_i +=  qp_data.ShapeValue(i, qp) * qp_data.JxW(qp);
      vol_frac_shape_i /= cell_mapping.CellVolume();
      VF[counter][i] = vol_frac_shape_i;
    }//for i

    counter++;
  }//for cell

}
