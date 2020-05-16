#include "chi_ffinter_line.h"
#include <ChiMath/SpatialDiscretization/PiecewiseLinear/pwl.h>

#include <chi_log.h>
extern ChiLog&  chi_log;

//###################################################################
/**Executes the interpolation.*/
void chi_mesh::FieldFunctionInterpolationLine::Execute()
{
  chi_log.Log(LOG_0VERBOSE_1) << "Executing line interpolator2.";
  for (int ff=0; ff<field_functions.size(); ff++)
  {
    grid_view = field_functions[ff]->grid;
    FieldFunctionContext* ff_ctx = ff_contexts[ff];

    if (field_functions[ff]->type == chi_physics::FieldFunctionType::CFEM_PWL)
    {
      Vec x_mapped;
      std::vector<int> mapping;
      Vec x = *field_functions[ff]->field_vector;
      CreateCFEMMapping(field_functions[ff]->num_components,
                        field_functions[ff]->num_sets,
                        field_functions[ff]->ref_component,
                        field_functions[ff]->ref_set,
                        x,x_mapped,
                        ff_ctx->cfem_local_nodes_needed_unmapped,
                        &mapping,
                        field_functions[ff]->spatial_discretization);

      CFEMInterpolate(x_mapped,mapping,ff_ctx);

    }
    else if (field_functions[ff]->type == chi_physics::FieldFunctionType::DFEM_PWL)
    {
      std::vector<int> mapping;
      CreatePWLDMapping(field_functions[ff]->num_components,
                        field_functions[ff]->num_sets,
                        field_functions[ff]->ref_component,
                        field_functions[ff]->ref_set,
                        ff_ctx->pwld_local_nodes_needed_unmapped,
                        ff_ctx->pwld_local_cells_needed_unmapped,
                        field_functions[ff]->spatial_discretization->cell_dfem_block_address,
                        &mapping);
      PWLDInterpolate(mapping,ff_ctx);
    }
    else if (field_functions[ff]->type == chi_physics::FieldFunctionType::FV)
    {
      std::vector<int> mapping;
      CreateFVMapping(field_functions[ff]->num_components,
                      field_functions[ff]->num_sets,
                      field_functions[ff]->ref_component,
                      field_functions[ff]->ref_set,
                      ff_ctx->interpolation_points_ass_cell,
                      &mapping);
      FVInterpolate(mapping,ff_ctx);
    }
  }

}

//###################################################################
/**Computes the cell average of each cell that was cut.*/
void chi_mesh::FieldFunctionInterpolationLine::
  FVInterpolate(std::vector<int> &mapping,
                FieldFunctionContext *ff_ctx)
{
  ff_ctx->interpolation_points_values.resize(interpolation_points.size(),0.0);

  for (size_t c=0; c<ff_ctx->interpolation_points_ass_cell.size(); c++)
  {
    if (ff_ctx->interpolation_points_ass_cell[c] < 0) continue;

    ff_ctx->interpolation_points_values[c] =
      ff_ctx->ref_ff->field_vector_local->operator[](mapping[c]);
  }
}

//###################################################################
/**Computes the cell average of each cell that was cut.*/
void chi_mesh::FieldFunctionInterpolationLine::
CFEMInterpolate(Vec field,
                std::vector<int> &mapping,
                FieldFunctionContext* ff_ctx)
{
  auto spatial_dm   =
    (SpatialDiscretization_PWL*)ff_ctx->ref_ff->spatial_discretization;

  //================================================== Loop over node indices
  //                                                   that
  //                                                   need mapping
  ff_ctx->interpolation_points_values.resize(interpolation_points.size(),0.0);
  int counter = -1;
  for (int c=0; c<ff_ctx->interpolation_points_ass_cell.size(); c++)
  {
    if (ff_ctx->interpolation_points_ass_cell[c] < 0) continue;

    int cell_local_index = ff_ctx->interpolation_points_ass_cell[c];
    auto cell_fe_view = spatial_dm->MapFeViewL(cell_local_index);

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
  PWLDInterpolate(std::vector<int> &mapping,
                  FieldFunctionContext* ff_ctx)
{
  auto spatial_dm   =
    (SpatialDiscretization_PWL*)ff_ctx->ref_ff->spatial_discretization;

  std::vector<double>& field = *ff_ctx->ref_ff->field_vector_local;

  //================================================== Loop over node indices
  //                                                   that
  //                                                   need mapping
  ff_ctx->interpolation_points_values.resize(interpolation_points.size(),0.0);
  int counter = -1;
  for (int c=0; c<ff_ctx->interpolation_points_ass_cell.size(); c++)
  {
    if (ff_ctx->interpolation_points_ass_cell[c] < 0) continue;

    int cell_local_index = ff_ctx->interpolation_points_ass_cell[c];
    auto cell_fe_view = spatial_dm->MapFeViewL(cell_local_index);

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