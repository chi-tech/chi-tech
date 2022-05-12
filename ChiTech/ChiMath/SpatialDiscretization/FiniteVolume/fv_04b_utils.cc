#include "fv.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>

#include "ChiMath/UnknownManager/unknown_manager.h"
#include "ChiMath/PETScUtils/petsc_utils.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;


//###################################################################
/**Get the number of local degrees-of-freedom.*/
size_t SpatialDiscretization_FV::
  GetNumLocalDOFs(chi_math::UnknownManager& unknown_manager)
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return ref_grid->local_cells.size()*N;
}

//###################################################################
/**Get the number of global degrees-of-freedom.*/
size_t SpatialDiscretization_FV::
  GetNumGlobalDOFs(chi_math::UnknownManager& unknown_manager)
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  const int num_globl_cells = ref_grid->GetGlobalNumberOfCells();

  return num_globl_cells*N;
}

//###################################################################
/**Get the number of ghost degrees-of-freedom.*/
size_t SpatialDiscretization_FV::
  GetNumGhostDOFs(chi_mesh::MeshContinuumPtr grid,
                  chi_math::UnknownManager& unknown_manager)
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return grid->cells.GetNumGhosts()*N;
}

//###################################################################
/**Returns the ghost DOF indices.*/
std::vector<int> SpatialDiscretization_FV::
  GetGhostDOFIndices(chi_mesh::MeshContinuumPtr grid,
                     chi_math::UnknownManager& unknown_manager,
                     unsigned int unknown_id)
{
  std::vector<int> dof_ids;
  std::vector<uint64_t> ghost_cell_ids = grid->cells.GetGhostGlobalIDs();

  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  dof_ids.reserve(GetNumGhostDOFs(grid,unknown_manager));
  for (int cell_id : ghost_cell_ids)
    for (int comp=0; comp<N; ++comp)
    {
      dof_ids.push_back(
        MapDOF(grid->cells[cell_id],0,unknown_manager,unknown_id,comp) );
    }

  return dof_ids;
}

//###################################################################
/**Develops a localized view of a petsc vector.*/
void SpatialDiscretization_FV::
  LocalizePETScVector(Vec petsc_vector,
                      std::vector<double>& local_vector,
                      chi_math::UnknownManager& unknown_manager)
{
  size_t num_local_dofs = GetNumLocalDOFs(unknown_manager);

  chi_math::PETScUtils::CopyVecToSTLvector(petsc_vector,
                                           local_vector,
                                           num_local_dofs);
}