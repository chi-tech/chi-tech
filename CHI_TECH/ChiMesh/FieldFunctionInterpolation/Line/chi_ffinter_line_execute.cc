#include "chi_ffinter_line.h"
#include <ChiMath/SpatialDiscretization/PiecewiseLinear/pwl.h>

#include <chi_log.h>
extern ChiLog&  chi_log;

//###################################################################
/**Executes the interpolation.*/
void chi_mesh::FieldFunctionInterpolationLine::Execute()
{
  chi_log.Log(LOG_0VERBOSE_1) << "Executing line interpolator.";
  for (int ff=0; ff<field_functions.size(); ff++)
  {
    grid_view = field_functions[ff]->grid;
    FieldFunctionContext* ff_ctx = ff_contexts[ff];

    typedef chi_math::SpatialDiscretizationType SDMType;
    const auto& field_sdm_type = field_functions[ff]->spatial_discretization->type;

    if (field_sdm_type == SDMType::PIECEWISE_LINEAR_CONTINUOUS)
    {
      std::vector<std::pair<uint64_t,uint>> node_component_pairs;

      for (auto node_id : ff_ctx->cfem_local_nodes_needed_unmapped)
        node_component_pairs.emplace_back(node_id,ff_ctx->ref_ff->ref_component);

      Vec x_mapped;
      std::vector<uint64_t> mapping;

      ff_ctx->ref_ff->CreateCFEMMappingLocal(x_mapped,
                                             node_component_pairs,
                                             mapping);

      CFEMInterpolate(x_mapped,mapping,ff_ctx);

      VecDestroy(&x_mapped);

    }
    else if (field_sdm_type == SDMType::PIECEWISE_LINEAR_DISCONTINUOUS)
    {
      std::vector<std::tuple<uint64_t,uint,uint>> cell_node_component_tuples;

      size_t num_mappings = ff_ctx->pwld_local_cells_needed_unmapped.size();
      for (size_t m=0; m<num_mappings; ++m)
        cell_node_component_tuples.emplace_back(
          ff_ctx->pwld_local_cells_needed_unmapped[m],
          ff_ctx->pwld_local_nodes_needed_unmapped[m],
          ff_ctx->ref_ff->ref_component);

      std::vector<uint64_t> mapping;

      ff_ctx->ref_ff->CreatePWLDMappingLocal(cell_node_component_tuples,
                                             mapping);

      PWLDInterpolate(mapping,ff_ctx);
    }
    else if (field_sdm_type == SDMType::FINITE_VOLUME)
    {
      std::vector<std::pair<uint64_t,uint>> cell_component_pairs;

      for (int p=0; p<number_of_points; ++p)
        if (ff_ctx->interpolation_points_has_ass_cell[p])
          cell_component_pairs.emplace_back(
            ff_ctx->interpolation_points_ass_cell[p],
            ff_ctx->ref_ff->ref_component);

      std::vector<uint64_t> mapping;

      ff_ctx->ref_ff->CreateFVMappingLocal(cell_component_pairs,
                                           mapping);

      FVInterpolate(mapping,ff_ctx);
    }
  }

}

//###################################################################
/**Computes the cell average of each cell that was cut.*/
void chi_mesh::FieldFunctionInterpolationLine::
  FVInterpolate(std::vector<uint64_t>& mapping,
                FieldFunctionContext *ff_ctx)
{
  ff_ctx->interpolation_points_values.resize(interpolation_points.size(),0.0);

  int counter=-1;
  for (size_t c=0; c<ff_ctx->interpolation_points_ass_cell.size(); c++)
  {
    if (not ff_ctx->interpolation_points_has_ass_cell[c]) continue;

    ++counter;
    ff_ctx->interpolation_points_values[c] =
      ff_ctx->ref_ff->field_vector_local->operator[](mapping[counter]);
  }
}

//###################################################################
/**Computes the cell average of each cell that was cut.*/
void chi_mesh::FieldFunctionInterpolationLine::
  CFEMInterpolate(Vec field,
                  std::vector<uint64_t>& mapping,
                  FieldFunctionContext* ff_ctx)
{
  auto& spatial_dm = static_cast<SpatialDiscretization_PWL&>(
                     *ff_ctx->ref_ff->spatial_discretization);

  //================================================== Loop over node indices
  //                                                   that
  //                                                   need mapping
  ff_ctx->interpolation_points_values.resize(interpolation_points.size(),0.0);
  int counter = -1;
  for (int c=0; c<ff_ctx->interpolation_points_ass_cell.size(); c++)
  {
    if (not ff_ctx->interpolation_points_has_ass_cell[c]) continue;

    int cell_local_index = ff_ctx->interpolation_points_ass_cell[c];
    auto cell_fe_view = spatial_dm.MapFeViewL(cell_local_index);

    double weighted_value = 0.0;
    for (int i=0; i<cell_fe_view->dofs; i++)
    {
      double node_value=0.0;
      counter++;
      int ir = mapping[counter];
      VecGetValues(field,1,&ir,&node_value);

      double weight=0.0;
      //Here I use c in interpolation_points because the vector should
      //be one-to-one with it.
      weight = cell_fe_view->ShapeValue(i, interpolation_points[c]);

      node_value *= weight;

      weighted_value += node_value;
    }

    ff_ctx->interpolation_points_values[c] = weighted_value;
  }//for ass cell

}


//###################################################################
/**Computes the cell average of each cell that was cut.*/
void chi_mesh::FieldFunctionInterpolationLine::
  PWLDInterpolate(std::vector<uint64_t>& mapping,
                  FieldFunctionContext* ff_ctx)
{
  auto& spatial_dm = static_cast<SpatialDiscretization_PWL&>(
                     *ff_ctx->ref_ff->spatial_discretization);

  std::vector<double>& field = *ff_ctx->ref_ff->field_vector_local;

  //================================================== Loop over node indices
  //                                                   that
  //                                                   need mapping
  ff_ctx->interpolation_points_values.resize(interpolation_points.size(),0.0);
  int counter = -1;
  for (int c=0; c<ff_ctx->interpolation_points_ass_cell.size(); c++)
  {
    if (not ff_ctx->interpolation_points_has_ass_cell[c]) continue;

    int cell_local_index = ff_ctx->interpolation_points_ass_cell[c];
    auto cell_fe_view = spatial_dm.MapFeViewL(cell_local_index);

    double weighted_value = 0.0;
    for (int i=0; i<cell_fe_view->dofs; i++)
    {
      double node_value=0.0;
      counter++;
      int ir = mapping[counter];
      node_value = field[ir];

      double weight=0.0;
      //Here I use c in interpolation_points because the vector should
      //be one-to-one with it.
      weight = cell_fe_view->ShapeValue(i, interpolation_points[c]);

      node_value *= weight;

      weighted_value += node_value;
    }

    ff_ctx->interpolation_points_values[c] = weighted_value;
  }//for ass cell

}