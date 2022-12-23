#include "mg_diffusion_solver.h"
#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

//========================================================== Solve 1g problem
void mg_diffusion::Solver::Assemble_RHS(const unsigned int g)
{
  // copy the external source vector for group g into b
  VecSet(b, 0.0);
  VecCopy(bext[g], b);
  const double* xlocal;
  VecGetArrayRead(x[g], &xlocal);

  const auto& sdm  = *mg_diffusion::Solver::sdm_ptr;
  // compute inscattering term
  for (const auto& cell :  mg_diffusion::Solver::grid_ptr->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto  qp_data      = cell_mapping.MakeVolumeQuadraturePointData();
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto& xs = matid_to_xs_map.at(cell.material_id);
    const auto& S = xs->transfer_matrices[0];
    // const unsigned int ell = 0;

    for (size_t i=0; i<num_nodes; ++i)
    {
      const int64_t imap = sdm.MapDOF(cell,i);
      double inscatter_g = 0.0;

      for (size_t j = 0; j < num_nodes; ++j) {
        const int64_t jmap = sdm.MapDOFLocal(cell, j); //only serial
        for (const auto& [row_g, gprime, sigma_sm] : S.Row(g))
        {
          if (gprime != g) // jcr row_g ??
          {
            // get flux at node j
            const double flxj_gp = xlocal[jmap];
            for (size_t qp: qp_data.QuadraturePointIndices())
              inscatter_g += sigma_sm * flxj_gp *
                         qp_data.ShapeValue(i, qp) * qp_data.ShapeValue(j, qp) *
                         qp_data.JxW(qp);
          }//if gp!=g
        }// for gprime
      }//for j
      // add inscattering value to vector
      VecSetValue(b, imap, inscatter_g, ADD_VALUES);
    }//for i
  }//for cell
  VecRestoreArrayRead(x[g], &xlocal);

  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

}