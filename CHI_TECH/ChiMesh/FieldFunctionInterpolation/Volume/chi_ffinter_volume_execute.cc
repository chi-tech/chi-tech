#include "chi_ffinter_volume.h"
#include <ChiMesh/Cell/cell.h>
#include <ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h>

#include <chi_mpi.h>

//###################################################################
/**Executes the volume interpolation.*/
void chi_mesh::FieldFunctionInterpolationVolume::Execute()
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
void chi_mesh::FieldFunctionInterpolationVolume::
  CFEMInterpolate(Vec field, std::vector<uint64_t> &mapping)
{
  auto& discretization = static_cast<SpatialDiscretization_PWLD&>(
                         *field_functions[0]->spatial_discretization);

  int counter=-1;
  double total_volume = 0.0;
  op_value = 0.0;
  double max_value = 0.0;
  bool max_set = false;

  for (auto& cell : grid_view->local_cells)
  {
    bool inside_logvolume=true;

    if (logical_volume != nullptr)
      inside_logvolume = logical_volume->Inside(cell.centroid);

    if (inside_logvolume)
    {
      const auto& fe_intgrl_values = discretization.GetUnitIntegrals(cell);

      for (int i=0; i<cell.vertex_ids.size(); i++)
      {
        double value = 0.0;
        int ir = -1;

        counter++;
        ir = mapping[counter];
        VecGetValues(field,1,&ir,&value);

        if ((op_type >= OP_SUM_LUA) and (op_type <= OP_MAX_LUA))
          value = CallLuaFunction(value,cell.material_id);

        op_value += value*fe_intgrl_values.FIntV_shapeI(i);
        total_volume += fe_intgrl_values.FIntV_shapeI(i);

        if (!max_set)
        {
          max_value = value;
          max_set = true;
        }
        else
        {
          if (value > max_value)
            max_value = value;
        }
      }//for dof
    }//if inside logicalVol

  }//for local cell

  double all_value=0.0;
  double all_total_volume = 0.0;
  double all_max_value=0.0;

  MPI_Allreduce(&op_value,&all_value,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&total_volume,&all_total_volume,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&max_value,&all_max_value,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (op_type == OP_AVG)
    op_value = all_value/total_volume;

  if (op_type == OP_MAX)
    op_value = all_max_value;
}



//###################################################################
/**Computes the cell average of each cell that was cut.*/
void chi_mesh::FieldFunctionInterpolationVolume::
  PWLDInterpolate(std::vector<double>& field, std::vector<uint64_t> &mapping)
{
  auto& discretization = static_cast<SpatialDiscretization_PWLD&>(
                         *field_functions[0]->spatial_discretization);

  int counter=-1;
  op_value = 0.0;
  double max_value = 0.0;
  bool max_set = false;
  double total_volume = 0.0;

  for (auto& cell : grid_view->local_cells)
  {
    bool inside_logvolume=true;

    if (logical_volume != nullptr)
      inside_logvolume = logical_volume->Inside(cell.centroid);

    if (inside_logvolume)
    {
      const auto& fe_intgrl_values = discretization.GetUnitIntegrals(cell);

      for (int i=0; i < cell.vertex_ids.size(); i++)
      {
        double value = 0.0;
        int ir = -1;

        counter++;
        ir = mapping[counter];
        value = field[ir];

        if ((op_type >= OP_SUM_LUA) and (op_type <= OP_MAX_LUA))
          value = CallLuaFunction(value,cell.material_id);

        op_value += value*fe_intgrl_values.FIntV_shapeI(i);
        total_volume += fe_intgrl_values.FIntV_shapeI(i);

        if (!max_set)
        {
          max_value = value;
          max_set = true;
        }
        else
        {
          if (value > max_value)
            max_value = value;
        }
      }//for dof
    }//if inside logicalVol

  }//for local cell

  double all_value=0.0;
  double all_total_volume = 0.0;
  double all_max_value;

  MPI_Allreduce(&op_value,&all_value,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&total_volume,&all_total_volume,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&max_value,&all_max_value,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (op_type == OP_AVG)
    op_value = all_value/total_volume;
  else
    op_value = all_value;

  if (op_type == OP_MAX)
    op_value = all_max_value;
}