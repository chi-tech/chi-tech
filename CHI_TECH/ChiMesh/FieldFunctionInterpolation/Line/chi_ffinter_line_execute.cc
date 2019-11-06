#include "chi_ffinter_line.h"
#include "ChiMesh/Cell/cell_slabv2.h"
#include "ChiMesh/Cell/cell_polygonv2.h"
#include "ChiMesh/Cell/cell_polyhedronv2.h"
#include <ChiMath/SpatialDiscretization/PiecewiseLinear/pwl.h>
#include <ChiMath/SpatialDiscretization/PiecewiseLinear/CellViews/pwl_slab.h>
#include <ChiMath/SpatialDiscretization/PiecewiseLinear/CellViews/pwl_polygon.h>
#include <ChiMath/SpatialDiscretization/PiecewiseLinear/CellViews/pwl_polyhedron.h>

//###################################################################
/**Executes the interpolation.*/
void chi_mesh::FieldFunctionInterpolationLine::Execute()
{
  chi_log.Log(LOG_0VERBOSE_1) << "Executing line interpolator2.";
  for (int ff=0; ff<field_functions.size(); ff++)
  {
    grid_view = field_functions[ff]->grid;
    FieldFunctionContext* ff_ctx = ff_contexts[ff];

    if (field_functions[ff]->type == FF_SDM_CFEM)
    {
      Vec x_mapped;
      std::vector<int> mapping;
      Vec x = *field_functions[ff]->field_vector;
      CreateCFEMMapping(field_functions[ff]->num_grps,
                        field_functions[ff]->num_moms,
                        field_functions[ff]->grp,
                        field_functions[ff]->mom,
                        x,x_mapped,
                        ff_ctx->cfem_local_nodes_needed_unmapped,
                        &mapping);

      CFEMInterpolate(x_mapped,mapping,ff_ctx);

    }
    else if (field_functions[ff]->type == FF_SDM_PWLD)
    {
      std::vector<int> mapping;
      CreatePWLDMapping(field_functions[ff]->num_grps,
                        field_functions[ff]->num_moms,
                        field_functions[ff]->grp,
                        field_functions[ff]->mom,
                        ff_ctx->pwld_local_nodes_needed_unmapped,
                        ff_ctx->pwld_local_cells_needed_unmapped,
                        *field_functions[ff]->local_cell_dof_array_address,
                        &mapping);
      PWLDInterpolate(mapping,ff_ctx);
    }
    else if (field_functions[ff]->type == FF_SDM_FV)
    {
      std::vector<int> mapping;
      CreateFVMapping(field_functions[ff]->num_grps,
                      field_functions[ff]->num_moms,
                      field_functions[ff]->grp,
                      field_functions[ff]->mom,
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

  for (int c=0; c<ff_ctx->interpolation_points_ass_cell.size(); c++)
  {
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
  SpatialDiscretization_PWL* spatial_dm   =
    (SpatialDiscretization_PWL*)ff_ctx->ref_ff->spatial_discretization;

  //================================================== Loop over node indices
  //                                                   that
  //                                                   need mapping
  ff_ctx->interpolation_points_values.resize(interpolation_points.size(),0.0);
  int counter = -1;
  for (int c=0; c<ff_ctx->interpolation_points_ass_cell.size(); c++)
  {
    if (ff_ctx->interpolation_points_ass_cell[c] < 0) continue;

    int cell_glob_index = ff_ctx->interpolation_points_ass_cell[c];
    auto cell = grid_view->cells[cell_glob_index];

    if (cell->Type() == chi_mesh::CellType::CELL_NEWBASE)
    {
      auto cell_base = (chi_mesh::CellBase*)cell;

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
      if (cell_base->Type2() == chi_mesh::CellType::SLABV2)
      {
        auto slab_cell = (chi_mesh::CellSlabV2*)cell_base;
        SlabFEView*      cell_fe_view =
          (SlabFEView*)spatial_dm->MapFeView(cell_glob_index);

        double weighted_value = 0.0;
        for (int i=0; i<2; i++)
        {
          double node_value=0.0;
          counter++;
          int ir = mapping[counter];
          VecGetValues(field,1,&ir,&node_value);

          double weight=0.0;
          //Here I use c in interpolation_points because the vector should
          //be one-to-one with it.
          weight = cell_fe_view->Shape_x(i, interpolation_points[c]);

          node_value *= weight;

          weighted_value += node_value;
        }

        ff_ctx->interpolation_points_values[c] = weighted_value;
      }

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
      else if (cell_base->Type2() == chi_mesh::CellType::POLYGONV2)
      {
        auto poly_cell = (chi_mesh::CellPolygonV2*)cell_base;
        PolygonFEView*      cell_fe_view =
          (PolygonFEView*)spatial_dm->MapFeView(cell_glob_index);

        double weighted_value = 0.0;
        for (int i=0; i<poly_cell->vertex_ids.size(); i++)
        {
          double node_value=0.0;
          counter++;
          int ir = mapping[counter];
          VecGetValues(field,1,&ir,&node_value);

          double weight=0.0;
          //Here I use c in interpolation_points because the vector should
          //be one-to-one with it.
          weight = cell_fe_view->Shape_xy(i, interpolation_points[c]);

          node_value *= weight;

          weighted_value += node_value;
        }

        ff_ctx->interpolation_points_values[c] = weighted_value;
      }

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
      else if (cell_base->Type2() == chi_mesh::CellType::POLYHEDRONV2)
      {
        auto polyh_cell = (chi_mesh::CellPolyhedronV2*)cell_base;
        PolyhedronFEView*       cell_fe_view =
          (PolyhedronFEView*)spatial_dm->MapFeView(cell_glob_index);

        double weighted_value = 0.0;
        for (int i=0; i<polyh_cell->vertex_ids.size(); i++)
        {
          double node_value=0.0;
          counter++;
          int ir = mapping[counter];
          VecGetValues(field,1,&ir,&node_value);

          double weight=0.0;
          //Here I use c in interpolation_points because the vector should
          //be one-to-one with it.
          weight = cell_fe_view->Shape_xyz(i, interpolation_points[c]);

          node_value *= weight;

          weighted_value += node_value;
        }

        ff_ctx->interpolation_points_values[c] = weighted_value;

      }//if polyhedron
    }//new cell base

  }//for ass cell

}


//###################################################################
/**Computes the cell average of each cell that was cut.*/
void chi_mesh::FieldFunctionInterpolationLine::
  PWLDInterpolate(std::vector<int> &mapping,
                  FieldFunctionContext* ff_ctx)
{
  SpatialDiscretization_PWL* spatial_dm   =
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

    int cell_glob_index = ff_ctx->interpolation_points_ass_cell[c];
    auto cell = grid_view->cells[cell_glob_index];

    if (cell->Type() == chi_mesh::CellType::CELL_NEWBASE)
    {
      auto cell_base = (chi_mesh::CellBase*)cell;

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
      if (cell_base->Type2() == chi_mesh::CellType::SLABV2)
      {
        auto slab_cell = (chi_mesh::CellSlabV2*)cell_base;
        SlabFEView*      cell_fe_view =
          (SlabFEView*)spatial_dm->MapFeView(cell_glob_index);

        double weighted_value = 0.0;
        for (int i=0; i<2; i++)
        {
          double node_value=0.0;
          counter++;
          int ir = mapping[counter];
          node_value = field[ir];

          double weight=0.0;
          //Here I use c in interpolation_points because the vector should
          //be one-to-one with it.
          weight = cell_fe_view->Shape_x(i, interpolation_points[c]);

          node_value *= weight;

          weighted_value += node_value;
        }

        ff_ctx->interpolation_points_values[c] = weighted_value;
      }


      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
      if (cell_base->Type2() == chi_mesh::CellType::POLYHEDRONV2)
      {
        auto polyh_cell = (chi_mesh::CellPolyhedronV2*)cell_base;
        PolyhedronFEView*       cell_fe_view =
          (PolyhedronFEView*)spatial_dm->MapFeView(cell_glob_index);

        double weighted_value = 0.0;
        for (int i=0; i<polyh_cell->vertex_ids.size(); i++)
        {
          double node_value=0.0;
          counter++;
          int ir = mapping[counter];
          node_value = field[ir];

          double weight=0.0;
          //Here I use c in interpolation_points because the vector should
          //be one-to-one with it.
          weight = cell_fe_view->Shape_xyz(i, interpolation_points[c]);

          node_value *= weight;

          weighted_value += node_value;
        }

        ff_ctx->interpolation_points_values[c] = weighted_value;

      }//if polyh
    }//new cell base
  }//for ass cell

}