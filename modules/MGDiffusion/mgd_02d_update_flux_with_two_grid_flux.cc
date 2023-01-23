#include "mg_diffusion_solver.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

#include "chi_runtime.h"
#include "ChiLog/chi_log.h"
#include "ChiTimer/chi_timer.h"

//============================================= assemble matrix A
void mg_diffusion::Solver::Update_Flux_With_TwoGrid(const int64_t verbose)
{
  if (verbose > 2)
    chi::log.Log() << "\nUpdating Thermal fluxes from two-grid";

  const auto& grid = *grid_ptr;
  const auto& sdm  = *sdm_ptr;

  // contains two_grid flux, stored in last num_groups entry
  const double *xlocal_tg;
  VecGetArrayRead(x[num_groups], &xlocal_tg);

  int counter = 0;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& xstg = map_mat_id_2_tginfo.at(cell.material_id);

    for (unsigned int g = last_fast_group; g < num_groups; ++g)
    {
      for (size_t i=0; i<num_nodes; ++i)
      {
        const int64_t imap       = sdm.MapDOFLocal(cell, i);
        const int64_t imap_globl = sdm.MapDOF(cell, i);
        const double aux = xlocal_tg[imap] * VF[counter][i] * xstg.spectrum[g];
        VecSetValue(x[g], imap_globl, aux, ADD_VALUES);
      }// i
    }//g
    counter++;
  }//for cell

  // finalize
  for (unsigned int g = last_fast_group; g < num_groups; ++g)
  {
    VecAssemblyBegin(x[g]);
    VecAssemblyEnd(x[g]);
  }
  // release two-grid flux
  VecRestoreArrayRead(x[num_groups], &xlocal_tg);
}