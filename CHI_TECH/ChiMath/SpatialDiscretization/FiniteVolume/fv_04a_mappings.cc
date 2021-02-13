#include "fv.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Maps a finite volume degree of freedom. The default behavior is to
 * assume a nodal DOF storage scheme.*/
int SpatialDiscretization_FV::
MapDOF(chi_mesh::Cell* cell)
{
  if (cell == nullptr)
  {
    chi_log.Log(LOG_ALLERROR)
      << "SpatialDiscretization_FV::MapDOF reference cell is nullptr.";
    exit(EXIT_FAILURE);
  }

  size_t num_local_cells = ref_grid->local_cells.size();


  int address=-1;
  if (cell->partition_id == chi_mpi.location_id)
  {
    address = fv_local_block_address +
              cell->local_id;
  }
  else
  {
    int ghost_local_id = 0;
    for (auto ghost : neighbor_cells)
      if (ghost->global_id == cell->global_id)
        ghost_local_id = ghost->local_id;

    address = locJ_block_address[cell->partition_id] +
              ghost_local_id;
  }

  return address;
}

//###################################################################
/**Maps a finite volume degree of freedom using an unknown manager.*/
int SpatialDiscretization_FV::
MapDOF(chi_mesh::Cell* cell,
       chi_math::UnknownManager* unknown_manager,
       unsigned int unknown_id,
       unsigned int component)
{
  if (cell == nullptr)
  {
    chi_log.Log(LOG_ALLERROR)
      << "SpatialDiscretization_FV::MapDOF reference cell is nullptr.";
    exit(EXIT_FAILURE);
  }
  if (component < 0) return -1;

  auto storage = unknown_manager->dof_storage_type;

  size_t num_unknowns = unknown_manager->GetTotalUnknownStructureSize();
  size_t block_id     = unknown_manager->MapUnknown(unknown_id, component);
  size_t num_local_cells = ref_grid->local_cells.size();

  if (component >= num_unknowns) return -1;


  int address=-1;
  if (cell->partition_id == chi_mpi.location_id)
  {
    if (storage == chi_math::UnknownStorageType::BLOCK)
      address = fv_local_block_address*num_unknowns +
                num_local_cells*block_id +
                cell->local_id;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = fv_local_block_address*num_unknowns +
                cell->local_id*num_unknowns +
                block_id;
  }
  else
  {
    int ghost_local_id = 0;
    for (auto ghost : neighbor_cells)
      if (ghost->global_id == cell->global_id)
        ghost_local_id = ghost->local_id;

    if (storage == chi_math::UnknownStorageType::BLOCK)
      address = locJ_block_address[cell->partition_id] * num_unknowns +
                locJ_block_size[cell->partition_id] * block_id +
                ghost_local_id;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = locJ_block_address[cell->partition_id] * num_unknowns +
                ghost_local_id*num_unknowns +
                block_id;
  }

  return address;
}

//###################################################################
/**Maps a finite volume degree of freedom to a local address.
 * The default behavior is to assume a nodal DOF storage scheme.*/
int SpatialDiscretization_FV::
MapDOFLocal(chi_mesh::Cell* cell)
{
  if (cell == nullptr)
  {
    chi_log.Log(LOG_ALLERROR)
      << "SpatialDiscretization_FV::MapDOFLocal reference cell is nullptr.";
    exit(EXIT_FAILURE);
  }

  int address=-1;
  if (cell->partition_id == chi_mpi.location_id)
    address = cell->local_id;
  else
  {
    int ghost_local_id = 0;
    for (auto ghost : neighbor_cells)
      if (ghost->global_id == cell->global_id)
        ghost_local_id = ghost->local_id;

    address = ghost_local_id;
  }

  return address;
}

//###################################################################
/**Maps a finite volume degree of freedom to a local address using
 * an unknown manager.*/
int SpatialDiscretization_FV::
MapDOFLocal(chi_mesh::Cell* cell,
            chi_math::UnknownManager* unknown_manager,
            unsigned int unknown_id,
            unsigned int component)
{
  if (cell == nullptr)
  {
    chi_log.Log(LOG_ALLERROR)
      << "SpatialDiscretization_FV::MapDOFLocal reference cell is nullptr.";
    exit(EXIT_FAILURE);
  }
  if (component < 0) return -1;

  auto storage = unknown_manager->dof_storage_type;

  size_t num_unknowns = unknown_manager->GetTotalUnknownStructureSize();
  size_t block_id     = unknown_manager->MapUnknown(unknown_id, component);
  size_t num_local_cells = ref_grid->local_cells.size();

  if (component >= num_unknowns) return -1;


  int address=-1;
  if (cell->partition_id == chi_mpi.location_id)
  {
    if (storage == chi_math::UnknownStorageType::BLOCK)
      address = num_local_cells*block_id +
                cell->local_id;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = cell->local_id*num_unknowns +
                block_id;
  }
  else
  {
    int ghost_local_id = 0;
    for (auto ghost : neighbor_cells)
      if (ghost->global_id == cell->global_id)
        ghost_local_id = ghost->local_id;

    if (storage == chi_math::UnknownStorageType::BLOCK)
      address = locJ_block_size[cell->partition_id] * block_id +
                ghost_local_id;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = ghost_local_id*num_unknowns + block_id;

  }

  return address;
}