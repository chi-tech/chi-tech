#include "fieldfunction.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

#include <chi_log.h>
extern ChiLog& chi_log;

////###################################################################
///***/
//std::vector<double>& chi_physics::FieldFunction::
//  GetCellDOFValues(size_t cell_local_id, size_t component, size_t set)
//{
//  if (using_petsc_field_vector)
//  {
//    chi_log.Log(LOG_ALLERROR)
//      << "chi_physics::FieldFunction::GetCellDOFValues "
//      << "does not yet support petsc vectors.";
//    exit(EXIT_FAILURE);
//  }
//  else
//  {
//    auto cell = grid->local_cells[cell_local_id];
//    int num_dofs = cell.vertex_ids.size();
//
//    temp_cell_dof_values.resize(num_dofs,0.0);
//
//    int block_map = (*local_cell_dof_array_address)[cell_local_id];
//
//    for (int dof=0; dof<num_dofs; ++dof)
//    {
//      int ir = block_map + dof*num_components*num_sets +
//               num_components * set + component;
//      temp_cell_dof_values[dof] = (*field_vector_local)[ir];
//    }
//  }
//
//  return temp_cell_dof_values;
//}

//###################################################################
/**Computes mappings of cell-local id to unknown vector index.*/
void chi_physics::FieldFunction::
  CreateFVMappingLocal(std::vector<std::pair<uint64_t,uint>>& cell_component_pairs,
                       std::vector<uint64_t>& mapping)
{
  auto& sdm    = spatial_discretization;

  if (sdm->type != chi_math::SpatialDiscretizationType::FINITE_VOLUME)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Field function spatial"
                                " discretization is not of type "
                                " FINITE_VOLUME.");

  auto sdm_fv = std::static_pointer_cast<SpatialDiscretization_FV>(sdm);

  for (auto& cell_index_component_pair : cell_component_pairs)
  {
    uint64_t cell_local_index = cell_index_component_pair.first;
    unsigned int component    = cell_index_component_pair.second;

    auto& cell = grid->local_cells[cell_local_index];
    int address = sdm_fv->MapDOFLocal(cell,
                                      unknown_manager,
                                      ref_variable,
                                      component);

    mapping.push_back(address);
  }
}


//###################################################################
/** Creates a mapping from a global vector */
void chi_physics::FieldFunction::
CreateCFEMMappingLocal(Vec& x_mapped,
                       std::vector<std::pair<uint64_t,uint>>& node_component_pairs,
                       std::vector<uint64_t>& mapping)
{
  if (spatial_discretization->type !=
      chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_CONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Field function spatial"
                                " discretization is not of type "
                                " PIECEWISE_LINEAR_CONTINUOUS.");

  auto pwl_sdm = std::static_pointer_cast<SpatialDiscretization_PWLC>(spatial_discretization);

  size_t num_nodes_to_map = node_component_pairs.size();
  std::vector<int> mapped_nodes;
  for (size_t n=0; n< num_nodes_to_map; n++)
  {
    int ir = pwl_sdm->MapDOF(node_component_pairs[n].first,
                             unknown_manager,
                             ref_variable,
                             node_component_pairs[n].second);

    mapped_nodes.push_back(ir);
    mapping.push_back(n);
  }

  Vec x = *field_vector;

  VecCreateSeq(PETSC_COMM_SELF,mapped_nodes.size()+1,&x_mapped);
  VecSet(x_mapped, 0.0);

  std::vector<int> int_mapping;

  std::copy(mapping.begin(),
            mapping.end(),
            std::back_inserter(int_mapping));

  IS global_set;
  IS local_set;
  ISCreateGeneral(PETSC_COMM_WORLD, num_nodes_to_map, mapped_nodes.data(),
                  PETSC_COPY_VALUES,&global_set);
  ISCreateGeneral(PETSC_COMM_WORLD, num_nodes_to_map, int_mapping.data(),
                  PETSC_COPY_VALUES,&local_set);
  VecScatter scat;
  VecScatterCreate(x, global_set, x_mapped, local_set, &scat);
  VecScatterBegin(scat, x, x_mapped, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(scat, x, x_mapped, INSERT_VALUES, SCATTER_FORWARD);

  ISDestroy(&global_set);
  ISDestroy(&local_set);
}


//###################################################################
/**Computes interpolated field function values.*/
void chi_physics::FieldFunction::
CreatePWLDMappingLocal(
  std::vector<std::tuple<uint64_t,uint,uint>>& cell_node_component_tuples,
  std::vector<uint64_t>& mapping)
{
  if (spatial_discretization->type !=
      chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Field function spatial"
                                " discretization is not of type "
                                " PIECEWISE_LINEAR_DISCONTINUOUS.");

  auto pwl_sdm = std::static_pointer_cast<SpatialDiscretization_PWL>(spatial_discretization);

  for (const auto& data : cell_node_component_tuples)
  {
    uint64_t cell_local_index = std::get<0>(data);
    uint     node_number      = std::get<1>(data);
    uint     component_number = std::get<2>(data);

    auto& cell = grid->local_cells[cell_local_index];

    int address = pwl_sdm->MapDOFLocal(cell,
                                       node_number,
                                       unknown_manager,
                                       ref_variable,
                                       component_number);

    mapping.push_back(address);
  }//for each tuple

}