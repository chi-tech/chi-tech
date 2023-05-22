#include "chi_ffinter_slice.h"

#include "ChiMath/VectorGhostCommunicator/vector_ghost_communicator.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/finite_element.h"
#include "ChiPhysics/FieldFunction/fieldfunction_gridbased.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Executes the slice interpolation.*/
void chi_mesh::FieldFunctionInterpolationSlice::Execute()
{
  const auto& ref_ff = *field_functions_.front();
  const auto& sdm    = ref_ff.SDM();
  const auto& grid   = sdm.ref_grid_;

  const auto& uk_man = ref_ff.UnkManager();
  const auto uid = 0;
  const auto cid = ref_component_;

  const auto field_data = ref_ff.GetGhostedFieldVector();

  for (auto& cell_intersection : cell_intersections_)
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
