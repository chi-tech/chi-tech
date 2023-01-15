#include "mg_diffusion_solver.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

#include "chi_runtime.h"
#include "chi_log.h"

//========================================================== Solve 1g problem
void mg_diffusion::Solver::Assemble_RHS_TwoGrid(const int64_t verbose)
{
  if (verbose > 2)
    chi::log.Log() << "\nAssemblying RHS for two-grid ";

  VecSet(b, 0.0);

  const auto& sdm  = *sdm_ptr;
  // compute inscattering term
  for (const auto& cell :  grid_ptr->local_cells)
  {
    const auto &cell_mapping = sdm.GetCellMapping(cell);
    const auto qp_data = cell_mapping.MakeVolumeQuadraturePointData();
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto &S = matid_to_xs_map.at(cell.material_id)->transfer_matrices[0];

    for (unsigned g = last_fast_group; g < num_groups; ++g)
    {
      for (const auto &[row_g, gprime, sigma_sm]: S.Row(g)) {
        if (gprime > g) // the upper part for the residual of two-grid accel
        {
          const double *xlocal;
          const double *xlocal_old;
          VecGetArrayRead(x[gprime], &xlocal);
          VecGetArrayRead(x_old[gprime], &xlocal_old);

          for (size_t i = 0; i < num_nodes; ++i) {
            const int64_t imap = sdm.MapDOF(cell, i);
            double inscatter_g = 0.0;

            for (size_t j = 0; j < num_nodes; ++j) {
              const int64_t jmap = sdm.MapDOFLocal(cell, j);

              // get flux at node j
              const double delta_flxj_gp = xlocal[jmap] - xlocal_old[jmap];
              for (size_t qp: qp_data.QuadraturePointIndices())
                inscatter_g += sigma_sm * delta_flxj_gp *
                               qp_data.ShapeValue(i, qp) * qp_data.ShapeValue(j, qp) *
                               qp_data.JxW(qp);
            }//for j
            // add inscattering value to vector
            VecSetValue(b, imap, inscatter_g, ADD_VALUES);
          }//for i
          VecRestoreArrayRead(x[gprime]    , &xlocal);
          VecRestoreArrayRead(x_old[gprime], &xlocal_old);
        }//if gp!=g
      }// for gprime
    }// for g
  }//for cell

  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
//  VecView(b, PETSC_VIEWER_STDERR_WORLD);
//  chi::Exit(1234);

}
