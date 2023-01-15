#include "chi_ffinter_slice.h"

#include "ChiMath/VectorGhostCommunicator/vector_ghost_communicator.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"

//###################################################################
/**Executes the slice interpolation.*/
void chi_mesh::FieldFunctionInterpolationSlice::Execute()
{
  const auto& ref_ff = *field_functions.front();
  const auto& sdm    = *ref_ff.spatial_discretization;
  const auto& grid   = *sdm.ref_grid;

  const auto& uk_man = ref_ff.unknown_manager;
  const auto uid = ref_ff.ref_variable;
  const auto cid = ref_ff.ref_component;

  const auto field_data = ref_ff.GetGhostedFieldVector();

  for (auto& cell_intersection : cell_intersections)
  {
    const auto& cell = grid.local_cells[cell_intersection.ref_cell_local_id];
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    std::vector<double> dof_values(num_nodes, 0.0);
    for (size_t i=0; i<num_nodes; ++i)
    {
      const int64_t imap = sdm.MapDOFLocal(cell, i, uk_man, uid, cid);
      dof_values[i] = field_data[imap];
    }

    std::vector<double> shape_values(num_nodes, 0.0);
    for (auto& edge_intersection : cell_intersection.intersections)
    {
      cell_mapping.ShapeValues(edge_intersection.point, shape_values);
      double point_value = 0.0;
      for (size_t i=0; i<num_nodes; ++i)
        point_value += dof_values[i]*shape_values[i];

      edge_intersection.point_value = point_value;
    }//for edge intersection
  }//for cell intersection
}
