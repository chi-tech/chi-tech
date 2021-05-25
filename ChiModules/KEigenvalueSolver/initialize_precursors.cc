#include "k_eigenvalue_solver.h"

#include <ChiMesh/Cell/cell.h>
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

using namespace LinearBoltzmann;

//###################################################################
/**Computes the point wise delayed neutron precursor concentrations.*/
void KEigenvalue::Solver::InitializePrecursors()
{
  auto pwl = std::static_pointer_cast<SpatialDiscretization_PWLD>(discretization);

  Nj_new_local.assign(Nj_new_local.size(), 0.0);

  // ----- Loop over cells
  for (auto& cell : grid->local_cells)
  {
    // ----- Cell information
    const auto xs_id = matid_to_xs_map[cell.material_id];
    auto xs = material_xs[xs_id];
    auto cell_fe_view = pwl->GetCellMappingFE(cell.local_id);
    auto& transport_view = cell_transport_views[cell.local_id];

    // ----- Loop over cell dofs
    for (int i = 0; i < cell_fe_view->num_nodes; ++i)
    {
      int64_t ir = transport_view.MapDOF(i,0,0);
      int64_t jr = pwl->MapDOFLocal(cell, i, Nj_unk_man, 0, 0);
      double* Nj_newp  = &Nj_new_local[jr];
      double* phi_newp = &phi_new_local[ir];

      // ----- If a fissial material with precursors
      if ((xs->is_fissile) and (xs->num_precursors > 0)) {
        // ----- Loop over precursors
        for (int j = 0; j < xs->num_precursors; ++j)
        {
          int j_map = precursor_map[xs_id][j];

          // ----- Initialize precursors
          double coeff = xs->gamma[j] / xs->lambda[j];


          // ----- Loop over groups
          for (int g = 0; g < groups.size(); ++g)
            Nj_newp[j_map] += coeff * xs->nu_d_sigma_fg[g] *
                              phi_newp[g] / k_eff;
        }//for j
      }//if fissile and J > 0
    } //for i
  }//for c
}


