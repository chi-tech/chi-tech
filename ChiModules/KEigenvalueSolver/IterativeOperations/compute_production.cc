#include "../k_eigenvalue_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

using namespace LinearBoltzmann;

//###################################################################
/**Compute the total fission production in the problem.*/
double KEigenvalue::Solver::ComputeProduction()
{
  typedef SpatialDiscretization_FE  FE;
  const auto grid_fe_view = std::static_pointer_cast<FE>(discretization);

  int first_grp = groups.front().id;
  int last_grp = groups.back().id;

  //================================================== Loop over cells
  double local_production = 0.0;
  for (auto& cell : grid->local_cells)
  {
    auto& full_cell_view = cell_transport_views[cell.local_id];
    const auto& fe_intgrl_values = grid_fe_view->GetUnitIntegrals(cell);

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
      double* phi_newp = &phi_new_local[ir];
      double intV_shapeI = fe_intgrl_values.IntV_shapeI(i);

      // only fissile materials contribute
      if (xs->is_fissile)
      {
        //=================================== Loop over groups
        for (int g = first_grp; g <= last_grp; ++g)
        {
          //TODO: Once checks are in place to ensure nu = nu_prompt +
          //      nu_delayed, nu_sigma_f can be universally used here.
          double nu_sigma_f =
              (not options.use_precursors) ? xs->nu_sigma_f[g] :
              xs->nu_prompt_sigma_f[g] + xs->nu_delayed_sigma_f[g];

          local_production += nu_sigma_f * phi_newp[g] * intV_shapeI;
        }// for group
      }
    }//for node
  }//for cell

  double global_production = 0.0;
  MPI_Allreduce(&local_production, &global_production, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return global_production;
}