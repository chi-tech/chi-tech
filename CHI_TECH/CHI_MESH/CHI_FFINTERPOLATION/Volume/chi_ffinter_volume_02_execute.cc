#include "chi_ffinter_volume.h"
#include "../../CHI_CELL/cell_slab.h"
#include "../../CHI_CELL/cell_polygon.h"
#include "../../CHI_CELL/cell_polyhedron.h"
#include <CHI_MATH/CHI_DISCRETIZATION/CHI_DISCRETIZATION_PWL/pwl.h>
#include <CHI_MATH/CHI_DISCRETIZATION/CHI_DISCRETIZATION_PWL/CellViews/pwl_polyhedron.h>
#include <CHI_MATH/CHI_DISCRETIZATION/CHI_DISCRETIZATION_PWL/CellViews/pwl_polygon.h>
#include <CHI_MATH/CHI_DISCRETIZATION/CHI_DISCRETIZATION_PWL/CellViews/pwl_slab.h>

#include <chi_mpi.h>

//###################################################################
/**Executes the volume interpolation.*/
void chi_mesh::FieldFunctionInterpolationVolume::Execute()
{
  if (field_functions[0]->type == FF_SDM_CFEM)
  {
    Vec x_mapped;
    std::vector<int> mapping;
    Vec x = field_functions[0]->field_vector;
    CreateCFEMMapping(field_functions[0]->num_grps,
                      field_functions[0]->num_moms,
                      field_functions[0]->grp,
                      field_functions[0]->mom,
                      x,x_mapped,cfem_local_nodes_needed_unmapped,&mapping);

    CFEMInterpolate(x_mapped,mapping);

  }
  else if (field_functions[0]->type == FF_SDM_PWLD)
  {
    std::vector<int> mapping;
    CreatePWLDMapping(field_functions[0]->num_grps,
                      field_functions[0]->num_moms,
                      field_functions[0]->grp,
                      field_functions[0]->mom,
                      pwld_local_nodes_needed_unmapped,
                      pwld_local_cells_needed_unmapped,
                      *field_functions[0]->local_cell_dof_array_address,
                      &mapping);
    PWLDInterpolate(*field_functions[0]->field_vector_local,mapping);
  }
}

//###################################################################
/**Computes the cell average of each cell that was cut.*/
void chi_mesh::FieldFunctionInterpolationVolume::
CFEMInterpolate(Vec field, std::vector<int> &mapping)
{
  CHI_DISCRETIZATION_PWL* discretization =
    (CHI_DISCRETIZATION_PWL*) field_functions[0]->spatial_discretization;

  int counter=-1;
  double total_volume = 0.0;
  op_value = 0.0;
  size_t num_local_cells = grid_view->local_cell_glob_indices.size();
  for (int lc=0; lc<num_local_cells; lc++)
  {
    int cell_glob_index = grid_view->local_cell_glob_indices[lc];
    auto cell = grid_view->cells[cell_glob_index];

    bool inside_logvolume=true;

    if (logical_volume != nullptr)
      inside_logvolume = logical_volume->Inside(cell->centroid);

    if (inside_logvolume)
    {
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
      if (typeid(*cell) == typeid(chi_mesh::CellSlab))
      {
        chi_mesh::CellSlab* slab_cell = (chi_mesh::CellSlab*)cell;
        SlabFEView* cell_fe_view =
          (SlabFEView*)discretization->MapFeView(cell_glob_index);

        for (int i=0; i<2; i++)
        {
          double value = 0.0;
          int ir = -1;

          counter++;
          ir = mapping[counter];
          VecGetValues(field,1,&ir,&value);

          op_value += value*cell_fe_view->IntV_shapeI[i];
          total_volume += cell_fe_view->IntV_shapeI[i];


        }//for dof
      }//if slab

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
      if (typeid(*cell) == typeid(chi_mesh::CellPolygon))
      {
        chi_mesh::CellPolygon* poly_cell = (chi_mesh::CellPolygon*)cell;
        PolygonFEView* cell_fe_view =
          (PolygonFEView*)discretization->MapFeView(cell_glob_index);

        for (int i=0; i<poly_cell->v_indices.size(); i++)
        {
          double value = 0.0;
          int ir = -1;

          counter++;
          ir = mapping[counter];
          VecGetValues(field,1,&ir,&value);

          op_value += value*cell_fe_view->IntV_shapeI[i];
          total_volume += cell_fe_view->IntV_shapeI[i];


        }//for dof
      }//if polygon

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
      if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
      {
        chi_mesh::CellPolyhedron* polyh_cell = (chi_mesh::CellPolyhedron*)cell;
        PolyhedronFEView* cell_fe_view =
          (PolyhedronFEView*)discretization->MapFeView(cell_glob_index);

        for (int i=0; i<polyh_cell->v_indices.size(); i++)
        {
          double value = 0.0;
          int ir = -1;

          counter++;
          ir = mapping[counter];
          VecGetValues(field,1,&ir,&value);

          op_value += value*cell_fe_view->IntV_shapeI[i];
          total_volume += cell_fe_view->IntV_shapeI[i];


        }//for dof
      }//if Polyhedron
    }//if inside logicalVol

  }//for local cell

  double all_value=0.0;
  double all_total_volume = 0.0;

  MPI_Allreduce(&op_value,&all_value,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&total_volume,&all_total_volume,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if (op_type == OP_AVG)
    op_value = all_value/total_volume;
}



//###################################################################
/**Computes the cell average of each cell that was cut.*/
void chi_mesh::FieldFunctionInterpolationVolume::
PWLDInterpolate(std::vector<double>& field, std::vector<int> &mapping)
{
  CHI_DISCRETIZATION_PWL* discretization =
    (CHI_DISCRETIZATION_PWL*) field_functions[0]->spatial_discretization;

  int counter=-1;
  op_value = 0.0;
  double total_volume = 0.0;
  size_t num_local_cells = grid_view->local_cell_glob_indices.size();
  for (int lc=0; lc<num_local_cells; lc++)
  {
    int cell_glob_index = grid_view->local_cell_glob_indices[lc];
    auto cell = grid_view->cells[cell_glob_index];

    bool inside_logvolume=true;

    if (logical_volume != nullptr)
      inside_logvolume = logical_volume->Inside(cell->centroid);

    if (inside_logvolume)
    {
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
      if (typeid(*cell) == typeid(chi_mesh::CellPolyhedron))
      {
        chi_mesh::CellPolyhedron* polyh_cell = (chi_mesh::CellPolyhedron*)cell;
        PolyhedronFEView* cell_fe_view =
          (PolyhedronFEView*)discretization->MapFeView(cell_glob_index);

        for (int i=0; i<polyh_cell->v_indices.size(); i++)
        {
          double value = 0.0;
          int ir = -1;

          counter++;
          ir = mapping[counter];
          value = field[ir];

          op_value += value*cell_fe_view->IntV_shapeI[i];
          total_volume += cell_fe_view->IntV_shapeI[i];

        }//for dof
      }//if Polyhedron
    }//if inside logicalVol

  }//for local cell

  double all_value=0.0;
  double all_total_volume = 0.0;

  MPI_Allreduce(&op_value,&all_value,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&total_volume,&all_total_volume,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if (op_type == OP_AVG)
    op_value = all_value/total_volume;
  else
    op_value = all_value;
}