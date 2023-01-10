#include "lbsadj_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"


/**Main execution function.*/
void lbs_adjoint::AdjointSolver::Execute()
{
  const std::string fname = __FUNCTION__;
  lbs::SteadySolver::Execute();

  //============================================= Apply post processing
  chi::log.Log() << "LBAdjointSolver: post-processing.";
  std::set<int> set_group_numbers;
  for (const auto& groupset : groupsets)
    for (const auto& group : groupset.groups)
      set_group_numbers.insert(group.id);

  const auto& m_to_ell_em_map =
    groupsets.front().quadrature->GetMomentToHarmonicsIndexMap();

  //============================================= Reorient phi-moments for reverse
  //                                              angle
  for (const auto& cell : grid->local_cells)
  {
    const auto& cell_view = cell_transport_views[cell.local_id];
    const int num_nodes = cell_view.NumNodes();

    for (int i = 0; i < num_nodes; ++i)
    {
      for (int m = 0; m < num_moments; ++m)
      {
        const auto& ell = m_to_ell_em_map[m].ell;

        size_t dof_map_g0 = cell_view.MapDOF(i, m, 0); //unknown map

        for (int g : set_group_numbers)
          phi_old_local[dof_map_g0 + g] *= pow(-1.0,ell);
      }//for moment
    }//node i
  }//for cell

}