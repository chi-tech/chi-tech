#include "lbs_linear_boltzmann_solver.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"

//###################################################################
/**Copy relevant section of phi_old to the field functions.*/
void lbs::SteadySolver::UpdateFieldFunctions()
{
  const auto& sdm = *discretization;
  const auto& phi_uk_man = flux_moments_uk_man;

  int ff_index = 0;
  for (size_t g=0; g<groups.size(); ++g)
  {
    for (size_t m=0; m<num_moments; ++m)
    {
      std::vector<double> data_vector_local(local_node_count, 0.0);

      for (const auto& cell : grid->local_cells)
      {
        const auto& cell_mapping = sdm.GetCellMapping(cell);
        const size_t num_nodes = cell_mapping.NumNodes();

        for (size_t i=0; i<num_nodes; ++i)
        {
          const int64_t imapA = sdm.MapDOFLocal(cell,i,phi_uk_man,m,g);
          const int64_t imapB = sdm.MapDOFLocal(cell,i);

          data_vector_local[imapB] = phi_old_local[imapA];
        }//for node
      }//for cell

      auto& ff_ptr = field_functions.at(ff_index);
      ff_ptr->UpdateFieldVector(data_vector_local);
      ++ff_index;
    }//for m
  }//for g

}

