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
  if (options.use_precursors)
  {
    typedef SpatialDiscretization_PWLD  PWLD;
    auto pwl = std::static_pointer_cast<PWLD>(discretization);

    precursor_new_local.assign(precursor_new_local.size(), 0.0);

    //======================================== Loop over cells
    for (auto& cell : grid->local_cells)
    {
      //==================== Cell information
      const auto xs_id = matid_to_xs_map[cell.material_id];
      auto& xs = material_xs[xs_id];
      auto cell_fe_view = pwl->GetCellMappingFE(cell.local_id);
      auto& transport_view = cell_transport_views[cell.local_id];

      //============================== Loop over nodes
      for (int i = 0; i < cell_fe_view->num_nodes; ++i)
      {
        size_t ir = transport_view.MapDOF(i, 0, 0);
        size_t jr = pwl->MapDOFLocal(cell, i, precursor_uk_man, 0, 0);
        double* Nj_newp = &precursor_new_local[jr];
        double* phi_newp = &phi_new_local[ir];

        // Contribute if precursors live on this material
        if (xs->num_precursors > 0)
        {
          //======================================== Loop over precursors
          for (size_t j = 0; j < xs->num_precursors; ++j)
          {

            //======================================== Loop over groups
            for (int g = 0; g < groups.size(); ++g)
              Nj_newp[j] += xs->precursor_yield[j] /
                            xs->precursor_lambda[j] *
                            xs->nu_delayed_sigma_f[g] *
                            phi_newp[g] / k_eff;
          }//for precursors
        }//if num_precursors > 0
      } //for node
    }//for cell
  }
}


