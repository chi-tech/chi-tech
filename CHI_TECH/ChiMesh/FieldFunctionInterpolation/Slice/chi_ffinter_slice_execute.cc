#include "chi_ffinter_slice.h"

//###################################################################
/**Executes the slice interpolation.*/
void chi_mesh::FieldFunctionInterpolationSlice::Execute()
{
  auto& ref_ff = *field_functions.back();
  const auto& field_sdm_type = ref_ff.spatial_discretization->type;

  typedef chi_math::SpatialDiscretizationType SMDType;

  if (field_sdm_type == SMDType::PIECEWISE_LINEAR_CONTINUOUS)
  {
    std::vector<std::pair<uint64_t,uint>> node_component_pairs;

    for (auto node_id : cfem_local_nodes_needed_unmapped)
      node_component_pairs.emplace_back(node_id, ref_ff.ref_component);

    Vec x_mapped;
    std::vector<uint64_t> mapping;

    ref_ff.CreateCFEMMappingLocal(x_mapped,
                                  node_component_pairs,
                                  mapping);

    CFEMInterpolate(x_mapped,mapping);

  }
  else if (field_sdm_type == SMDType::PIECEWISE_LINEAR_DISCONTINUOUS)
  {
    std::vector<std::tuple<uint64_t,uint,uint>> cell_node_component_tuples;

    size_t num_mappings = pwld_local_cells_needed_unmapped.size();
    for (size_t m=0; m<num_mappings; ++m)
      cell_node_component_tuples.emplace_back(
        pwld_local_cells_needed_unmapped[m],
        pwld_local_nodes_needed_unmapped[m],
        ref_ff.ref_component);

    std::vector<uint64_t> mapping;

    ref_ff.CreatePWLDMappingLocal(cell_node_component_tuples,
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
