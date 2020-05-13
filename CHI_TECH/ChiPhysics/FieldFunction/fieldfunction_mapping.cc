#include "fieldfunction.h"

#include <chi_log.h>
extern ChiLog& chi_log;

//###################################################################
/***/
std::vector<double>& chi_physics::FieldFunction::
  GetCellDOFValues(size_t cell_local_id, size_t component, size_t set)
{
  if (using_petsc_field_vector)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_physics::FieldFunction::GetCellDOFValues "
      << "does not yet support petsc vectors.";
    exit(EXIT_FAILURE);
  }
  else
  {
    auto cell = grid->local_cells[cell_local_id];
    int num_dofs = cell.vertex_ids.size();

    temp_cell_dof_values.resize(num_dofs,0.0);

    int block_map = (*local_cell_dof_array_address)[cell_local_id];

    for (int dof=0; dof<num_dofs; ++dof)
    {
      int ir = block_map + dof*num_components*num_sets +
               num_components * set + component;
      temp_cell_dof_values[dof] = (*field_vector_local)[ir];
    }
  }

  return temp_cell_dof_values;
}