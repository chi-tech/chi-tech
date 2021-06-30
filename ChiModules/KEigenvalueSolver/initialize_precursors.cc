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
    precursor_new_local.assign(precursor_new_local.size(), 0.0);

    //================================================== Loop over cells
    for (auto& cell : grid->local_cells)
    {
      auto& full_cell_view = cell_transport_views[cell.local_id];

      //==================== Obtain xs
      int cell_matid = cell.material_id;
      int xs_id = matid_to_xs_map[cell_matid];

      if ((xs_id < 0) || (xs_id >= material_xs.size()))
      {
        chi_log.Log(LOG_ALLERROR)
            << "Cross-section lookup error\n";
        exit(EXIT_FAILURE);
      }

      auto xs = material_xs[xs_id];

      //======================================== Loop over nodes
      const int num_nodes = full_cell_view.NumNodes();
      for (int i = 0; i < num_nodes; ++i)
      {
        size_t ir = full_cell_view.MapDOF(i, 0, 0);
        size_t jr = pwl->MapDOFLocal(cell, i, precursor_uk_man, 0, 0);
        double* Nj_newp = &precursor_new_local[jr];
        double* phi_newp = &phi_new_local[ir];

        // contribute if precursors live on this material
        if (xs->num_precursors > 0)
        {
          //======================================== Loop over precursors
          for (size_t j = 0; j < xs->num_precursors; ++j)
          {

            //======================================== Loop over groups
            for (size_t g = 0; g < groups.size(); ++g)
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


