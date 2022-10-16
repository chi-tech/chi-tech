#include "fv.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell.h"

#include "chi_log.h"
#include "chi_mpi.h"

#define sc_int64 static_cast<int64_t>

//###################################################################
/**Maps a finite volume degree of freedom using an unknown manager.*/
int64_t chi_math::SpatialDiscretization_FV::
  MapDOF(const chi_mesh::Cell& cell,
         const unsigned int,
         const chi_math::UnknownManager& unknown_manager,
         const unsigned int unknown_id,
         const unsigned int component) const
{
  auto storage = unknown_manager.dof_storage_type;

  size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  size_t block_id     = unknown_manager.MapUnknown(unknown_id, component);
  size_t num_local_cells = ref_grid->local_cells.size();

  if (component >= num_unknowns) return -1;


  int64_t address=-1;
  if (cell.partition_id == chi::mpi.location_id)
  {
    if (storage == chi_math::UnknownStorageType::BLOCK)
      address = sc_int64(local_block_address) * num_unknowns +
                num_local_cells*block_id +
                cell.local_id;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = sc_int64(local_block_address) * num_unknowns +
                cell.local_id*num_unknowns +
                block_id;
  }
  else
  {
    uint64_t ghost_local_id = neighbor_cells.at(cell.global_id)->local_id;

    if (storage == chi_math::UnknownStorageType::BLOCK)
      address = sc_int64(locJ_block_address[cell.partition_id]) * num_unknowns +
                locJ_block_size[cell.partition_id] * block_id +
                ghost_local_id;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = sc_int64(locJ_block_address[cell.partition_id]) * num_unknowns +
                ghost_local_id*num_unknowns +
                block_id;
  }

  return address;
}

//###################################################################
/**Maps a finite volume degree of freedom to a local address using
 * an unknown manager.*/
int64_t chi_math::SpatialDiscretization_FV::
  MapDOFLocal(const chi_mesh::Cell& cell,
              const unsigned int,
              const chi_math::UnknownManager& unknown_manager,
              const unsigned int unknown_id,
              const unsigned int component) const
{
  auto storage = unknown_manager.dof_storage_type;

  size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  size_t block_id     = unknown_manager.MapUnknown(unknown_id, component);
  size_t num_local_cells = ref_grid->local_cells.size();

  if (component >= num_unknowns) return -1;


  int address=-1;
  if (cell.partition_id == chi::mpi.location_id)
  {
    if (storage == chi_math::UnknownStorageType::BLOCK)
      address = sc_int64(num_local_cells) * block_id + cell.local_id;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = sc_int64(cell.local_id) * num_unknowns + block_id;
  }
  else
  {
    uint64_t ghost_local_id = neighbor_cells.at(cell.global_id)->local_id;

    if (storage == chi_math::UnknownStorageType::BLOCK)
      address = sc_int64(locJ_block_size[cell.partition_id]) * block_id +
                ghost_local_id;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = sc_int64(ghost_local_id)*num_unknowns + block_id;

  }

  return address;
}