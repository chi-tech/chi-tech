#include "mg_diffusion_solver.h"
#include "ChiTimer/chi_timer.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

//============================================= assemble matrix A
void mg_diffusion::Solver::Update_Flux_With_TwoGrid(const int64_t verbose)
{
  const auto& grid = *grid_ptr;
  const auto& sdm  = *sdm_ptr;

  const double *xlocal_tg;
  VecGetArrayRead(x[num_groups], &xlocal_tg);

  int counter = 0;
  for (const auto& cell : grid.local_cells)
  {
    const double *xlocal;

    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& xstg = map_mat_id_2_tginfo.at(cell.material_id);

    for (unsigned int g = last_fast_group; g < num_groups; ++g)
    {
      VecGetArrayRead(x[g], &xlocal);
      for (size_t i=0; i<num_nodes; ++i)
      {
        const int64_t imap = sdm.MapDOFLocal(cell, i); //MapDOFLocal ??
        xlocal[imap] += xlocal_tg[imap] * VF[counter][i] * xstg.spectrum[g];
      }// i
      VecRestoreArrayRead(x[g], &xlocal);
    }//g
    counter++;
  }//for cell
  
  VecRestoreArrayRead(x[num_groups], &xlocal_tg);
}