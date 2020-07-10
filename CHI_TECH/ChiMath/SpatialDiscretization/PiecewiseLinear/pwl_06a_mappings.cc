#include "pwl.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Maps a vertex id according to a developed node ordering.*/
int SpatialDiscretization_PWL::MapCFEMDOF(int vertex_id)
{
  if (node_mapping.empty())
    return vertex_id;

  return node_mapping[vertex_id];
}

//###################################################################
/**Maps a vertex id according to a developed node ordering.*/
int SpatialDiscretization_PWL::
  MapCFEMDOF(int vertex_id,
             chi_math::UnknownManager* unknown_manager,
             unsigned int unknown_id,
             unsigned int component)
{
  if (component < 0) return -1;

  auto storage = unknown_manager->dof_storage_type;

  int mapping = vertex_id;
  if (not node_mapping.empty()) mapping = node_mapping[vertex_id];

  size_t num_unknowns = unknown_manager->GetTotalUnknownSize();
  size_t block_id     = unknown_manager->MapUnknown(unknown_id,component);

  int address=-1;
  if (storage == chi_math::DOFStorageType::BLOCK)
  {
    for (int locJ=0; locJ<chi_mpi.process_count; ++locJ)
    {
      int localized_base_address = mapping - locJ_block_address[locJ];
      if (localized_base_address < 0) continue;
      if (localized_base_address >= locJ_block_size[locJ]) continue;

      address = locJ_block_address[locJ]*num_unknowns +
                locJ_block_size[locJ]*block_id +
                localized_base_address;
      break;
    }
  }
  else if (storage == chi_math::DOFStorageType::NODAL)
    address = mapping*num_unknowns + block_id;

  return address;
}

//###################################################################
/**Provides a mapping of cell's DOF from a DFEM perspective.*/
int SpatialDiscretization_PWL::MapDFEMDOF(chi_mesh::Cell *cell,
                                          int dof,
                                          int component,
                                          int component_block_offset)
{
  if (cell->partition_id == chi_mpi.location_id)
  {
    int address = cell_dfem_block_address[cell->local_id] +
                  dfem_local_block_address +
                  dof;
    return address*component_block_offset + component;
  }
  else
  {
    int index = 0;
    bool found = false;
    for (auto neighbor_info : neighbor_cell_block_address)
    {
      if (neighbor_info.first == cell->global_id) {
        found = true; break;
      }
      ++index;
    }

    if (!found)
    {
      chi_log.Log(LOG_ALLERROR)
        << "SpatialDiscretization_PWL::MapDFEMDOF. Mapping failed for cell "
        << "with global index " << cell->global_id << " and partition-ID "
        << cell->partition_id;
      exit(EXIT_FAILURE);
    }

    int address = neighbor_cell_block_address[index].second + dof;
    return address*component_block_offset + component;
  }
}

//###################################################################
/**Provides a local mapping of cell's DOF from a DFEM perspective.*/
int SpatialDiscretization_PWL::MapDFEMDOFLocal(chi_mesh::Cell *cell,
                                               int dof,
                                               int component,
                                               int component_block_offset)
{
  if (cell->partition_id == chi_mpi.location_id)
  {
    int address = cell_dfem_block_address[cell->local_id] + dof;
    return address*component_block_offset + component;
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "SpatialDiscretization_PWL::MapDFEMDOF. Mapping failed for cell "
      << "with global index " << cell->global_id << " and partition-ID "
      << cell->partition_id;
    exit(EXIT_FAILURE);
  }
}

//###################################################################
/**Provides a mapping of cell's DOF from a DFEM perspective.*/
int SpatialDiscretization_PWL::
  MapDFEMDOF(chi_mesh::Cell *cell, int node,
             chi_math::UnknownManager* unknown_manager,
             unsigned int unknown_id,
             unsigned int component)
{
  if (component < 0) return -1;

  auto storage = unknown_manager->dof_storage_type;

  size_t num_unknowns = unknown_manager->GetTotalUnknownSize();
  size_t block_id     = unknown_manager->MapUnknown(unknown_id,component);

  if (cell->partition_id == chi_mpi.location_id)
  {
    if (storage == chi_math::DOFStorageType::BLOCK)
    {
      int address = dfem_local_block_address*num_unknowns +
                    cell_dfem_block_address[cell->local_id] +
                    local_base_block_size*block_id +
                    node;
      return address;
    }
    else if (storage == chi_math::DOFStorageType::NODAL)
    {
      int address = dfem_local_block_address*num_unknowns +
                    cell_dfem_block_address[cell->local_id]*num_unknowns +
                    node*num_unknowns +
                    block_id;
      return address;
    }
  }
  else
  {
    int index = 0;
    bool found = false;
    for (auto neighbor_info : neighbor_cell_block_address)
    {
      if (neighbor_info.first == cell->global_id) {
        found = true; break;
      }
      ++index;
    }

    if (!found)
    {
      chi_log.Log(LOG_ALLERROR)
        << "SpatialDiscretization_PWL::MapDFEMDOF. Mapping failed for cell "
        << "with global index " << cell->global_id << " and partition-ID "
        << cell->partition_id;
      exit(EXIT_FAILURE);
    }

    if (storage == chi_math::DOFStorageType::BLOCK)
    {
      int address = locJ_block_address[cell->partition_id]*num_unknowns +
                    neighbor_cell_block_address[index].second +
                    locJ_block_size[cell->partition_id]*block_id +
                    node;
      return address;
    }
    else if (storage == chi_math::DOFStorageType::NODAL)
    {
      int address = locJ_block_address[cell->partition_id]*num_unknowns +
                    neighbor_cell_block_address[index].second*num_unknowns +
                    node*num_unknowns +
                    block_id;
      return address;
    }

  }

  return -1;
}