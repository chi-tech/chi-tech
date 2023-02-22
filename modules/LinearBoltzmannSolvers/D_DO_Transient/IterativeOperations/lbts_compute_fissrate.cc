#include "D_DO_Transient/lbts_transient_solver.h"

using namespace lbs;

//###################################################################
/**Computes the total fission rate in the problem.

\param previous bool Optional. If true and the solver is a transient solver
                     then the rate will be based on the previous timestep's
                     fluxes.

\return the_rate The fission rate as a double.

\author Zachary Hardy.*/
double DiscOrdTransientSolver::ComputeFissionRate(const bool previous)
{
  const int first_grp = groups_.front().id;
  const int last_grp = groups_.back().id;

  const auto& phi = (previous) ? phi_prev_local_ : phi_new_local_;

  //============================================= Loop over local cells
  double local_production = 0.0;
  for (auto& cell : grid_ptr_->local_cells)
  {
    const auto& transport_view = cell_transport_views_[cell.local_id];
    const auto& cell_matrices = unit_cell_matrices_[cell.local_id];

    //====================================== Obtain xs
    const auto& xs = transport_view.XS();

    //====================================== Loop over nodes
    const int num_nodes = transport_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.Vi_vectors[i];

      //=============================== Loop over groups
      for (size_t g = first_grp; g <= last_grp; ++g)
        local_production += xs.sigma_f[g] *
                            phi[uk_map + g] *
                            IntV_ShapeI;
    }//for node
  }//for cell

  //============================================= Allreduce global production
  double global_production = 0.0;
  MPI_Allreduce(&local_production, &global_production, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return global_production;
}