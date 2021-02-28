#include "pwlc.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Maps a vertex id according to a developed node ordering.*/
int SpatialDiscretization_PWLC::
  MapDOF(const chi_mesh::Cell& cell,
         int node,
         const chi_math::UnknownManager& unknown_manager,
         unsigned int unknown_id,
         unsigned int component/*=0*/) const
{
  int vertex_id = cell.vertex_ids[node];

  auto storage = unknown_manager.dof_storage_type;

  int mapping = vertex_id;
  if (not node_mapping.empty()) mapping = node_mapping.at(vertex_id);

  size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  size_t block_id     = unknown_manager.MapUnknown(unknown_id, component);

  int address=-1;
  if (storage == chi_math::UnknownStorageType::BLOCK)
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
  else if (storage == chi_math::UnknownStorageType::NODAL)
    address = mapping*num_unknowns + block_id;

  return address;
}

