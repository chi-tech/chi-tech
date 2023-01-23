#include "chi_ffinter_line.h"

#include "ChiMath/VectorGhostCommunicator/vector_ghost_communicator.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"
#include "ChiMath/SpatialDiscretization/spatial_discretization.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Executes the interpolation.*/
void chi_mesh::FieldFunctionInterpolationLine::Execute()
{
  chi::log.Log0Verbose1() << "Executing line interpolator.";
  for (int ff=0; ff<field_functions.size(); ff++)
  {
          auto& ff_ctx = ff_contexts[ff];
    const auto& ref_ff = *ff_ctx.ref_ff;
    const auto& sdm    = ref_ff.SDM();
    const auto& grid   = *sdm.ref_grid;

    const auto& uk_man = ref_ff.UnkManager();
    const auto uid = 0;
    const auto cid = m_ref_component;

    const auto field_data = ref_ff.GetGhostedFieldVector();

    ff_ctx.interpolation_points_values.assign(number_of_points,0.0);
    for (int p=0; p<number_of_points; ++p)
    {
      if (not ff_ctx.interpolation_points_has_ass_cell[p]) continue;

      const auto cell_local_index = ff_ctx.interpolation_points_ass_cell[p];
      const auto& cell = grid.local_cells[cell_local_index];
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      std::vector<double> shape_function_vals(num_nodes, 0.0);
      cell_mapping.ShapeValues(interpolation_points[p], shape_function_vals);

      double point_value = 0.0;
      for (size_t i=0; i<num_nodes; ++i)
      {
        const int64_t imap = sdm.MapDOFLocal(cell, i, uk_man, uid, cid);

        point_value += shape_function_vals[i]*field_data[imap];
      }//for node i
      ff_ctx.interpolation_points_values[p] = point_value;
    }//for p
  }//for ff

}