#include "chi_ffinter_slice.h"

//###################################################################
/**Executes the slice interpolation.*/
void chi_mesh::FieldFunctionInterpolationSlice::Execute()
{
  if (field_functions[0]->type == chi_physics::FieldFunctionType::CFEM_PWL)
  {
    Vec x_mapped;
    std::vector<uint64_t> mapping;
    Vec x = *field_functions[0]->field_vector;

    CreateCFEMMapping(field_functions[0]->num_components,
                      field_functions[0]->num_sets,
                      field_functions[0]->ref_component,
                      field_functions[0]->ref_set,
                      x,x_mapped,cfem_local_nodes_needed_unmapped,mapping,
                      field_functions[0]->spatial_discretization);

    CFEMInterpolate(x_mapped,mapping);

  }
  else if (field_functions[0]->type == chi_physics::FieldFunctionType::DFEM_PWL)
  {
    std::vector<uint64_t> mapping;
    CreatePWLDMapping(field_functions[0]->num_components,
                      field_functions[0]->num_sets,
                      field_functions[0]->ref_component,
                      field_functions[0]->ref_set,
                      pwld_local_nodes_needed_unmapped,
                      pwld_local_cells_needed_unmapped,
                      field_functions[0]->spatial_discretization->cell_dfem_block_address,
                      mapping);
    PWLDInterpolate(*field_functions[0]->field_vector_local,mapping);
  }
}


//###################################################################
/**Computes the cell average of each cell that was cut.*/
void chi_mesh::FieldFunctionInterpolationSlice::
CFEMInterpolate(Vec field, std::vector<uint64_t> &mapping)
{
  size_t num_slice_cells = cell_intersections.size();
  int counter=-1;
  for (int sc=0; sc<num_slice_cells; sc++)
  {
    double cell_sum = 0.0;

    size_t num_is = cell_intersections[sc].intersections.size();
    for (int is=0; is<num_is; is++)
    {
      double value = 0.0;
      int ir = -1;

      counter++;
      ir = mapping[counter];
      VecGetValues(field,1,&ir,&value);

      cell_sum += value;

      cell_intersections[sc].intersections[is].point_value =
        cell_intersections[sc].intersections[is].weights.first*value;

      counter++;
      ir = mapping[counter];
      VecGetValues(field,1,&ir,&value);

      cell_sum += value;

      cell_intersections[sc].intersections[is].point_value +=
        cell_intersections[sc].intersections[is].weights.second*value;
    }

    cell_intersections[sc].cell_avg_value = cell_sum/2*num_is;
  }
}



//###################################################################
/**Computes the cell average of each cell that was cut.*/
void chi_mesh::FieldFunctionInterpolationSlice::
PWLDInterpolate(std::vector<double>& field, std::vector<uint64_t>& mapping)
{
  size_t num_slice_cells = cell_intersections.size();
  int counter=-1;
  for (int sc=0; sc<num_slice_cells; sc++)
  {
    double cell_sum = 0.0;

    size_t num_is = cell_intersections[sc].intersections.size();
    for (int is=0; is<num_is; is++)
    {
      double value = 0.0;
      int ir = -1;

      counter++;
      ir = mapping[counter];
      value = field[ir];


      cell_sum += value;

      cell_intersections[sc].intersections[is].point_value =
        cell_intersections[sc].intersections[is].weights.first*value;

      counter++;
      ir = mapping[counter];
      value = field[ir];

      cell_sum += value;

      cell_intersections[sc].intersections[is].point_value +=
        cell_intersections[sc].intersections[is].weights.second*value;
    }

    cell_intersections[sc].cell_avg_value = cell_sum/2*num_is;
  }
}
