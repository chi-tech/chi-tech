#include "lbs_linear_boltzmann_solver.h"

#include "ChiPhysics/FieldFunction2/fieldfunction2.h"

//###################################################################
/**Copy relevant section of phi_old to the field functions.*/
void lbs::SteadySolver::UpdateFieldFunctions()
{
  const auto& sdm = *discretization;
  const auto& phi_uk_man = flux_moments_uk_man;

  int    dimension = 0;
  if (grid->Attributes() & chi_mesh::DIMENSION_1) dimension = 1;
  if (grid->Attributes() & chi_mesh::DIMENSION_2) dimension = 2;
  if (grid->Attributes() & chi_mesh::DIMENSION_3) dimension = 3;

  std::array<unsigned int,4> m_map = {0, 0, 0, 0};
  if (dimension == 1 and num_moments >= 2) m_map = {0, 0, 0, 1};
  if (dimension == 2 and num_moments >= 3) m_map = {0, 2, 1, 0};
  if (dimension == 3 and num_moments >= 4) m_map = {0, 3, 1, 2};

  int ff_index = 0;
  for (size_t g=0; g<groups.size(); ++g)
  {
    for (size_t ff=0; ff<4; ++ff)
    {
      std::vector<double> data_vector_local(local_node_count, 0.0);

      const auto m = m_map[ff];
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

      auto& ff_ptr = field_functions2.at(ff_index);
      ff_ptr->UpdateFieldVector(data_vector_local);
      ++ff_index;
    }//for ff
  }//for g

}

