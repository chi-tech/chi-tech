#include "../lbkes_k_eigenvalue_solver.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

using namespace LinearBoltzmann;

//###################################################################
/**Compute the total fission production in the problem.*/
double KEigenvalueSolver::ComputeFissionProduction()
{
  typedef SpatialDiscretization_FE  FE;
  const auto grid_fe_view = std::dynamic_pointer_cast<FE>(discretization);

  int first_grp = groups.front().id;
  int last_grp = groups.back().id;

  //================================================== Loop over cells
  double local_production = 0.0;
  for (auto& cell : grid->local_cells)
  {
    const auto& transport_view = cell_transport_views[cell.local_id];
    const auto& fe_values = grid_fe_view->GetUnitIntegrals(cell);

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
    const int num_nodes = transport_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      size_t uk_map = transport_view.MapDOF(i, 0, 0);
      double IntV_ShapeI = fe_values.IntV_shapeI(i);

      //=================================== Loop over groups
      for (size_t g = first_grp; g <= last_grp; ++g)
        local_production += xs->nu_sigma_f[g] *
                            phi_new_local[uk_map + g] *
                            IntV_ShapeI;
    }//for node
  }//for cell

  double global_production = 0.0;
  MPI_Allreduce(&local_production, &global_production, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return global_production;
}