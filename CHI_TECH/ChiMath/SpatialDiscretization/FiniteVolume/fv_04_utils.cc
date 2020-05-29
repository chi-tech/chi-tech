#include "fv.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>

#include "ChiMath/UnknownManager/unknown_manager.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Get the number of local degrees-of-freedom.*/
unsigned int SpatialDiscretization_FV::
  GetNumLocalDOFs(chi_mesh::MeshContinuum* grid,
                  chi_math::UnknownManager* unknown_manager)
{
  unsigned int N = 1;

  if (unknown_manager != nullptr)
    N = unknown_manager->GetTotalUnknownSize();

  const int num_local_cells = grid->local_cells.size();

  return num_local_cells*N;
}

//###################################################################
/**Get the number of local degrees-of-freedom.*/
unsigned int SpatialDiscretization_FV::
  GetNumGlobalDOFs(chi_mesh::MeshContinuum* grid,
                   chi_math::UnknownManager* unknown_manager)
{
  unsigned int N = 1;

  if (unknown_manager != nullptr)
    N = unknown_manager->GetTotalUnknownSize();

  const int num_globl_cells = grid->GetGlobalNumberOfCells();

  return num_globl_cells*N;
}

//###################################################################
/**Get the number of ghost degrees-of-freedom.*/
unsigned int SpatialDiscretization_FV::
  GetNumGhostDOFs(chi_mesh::MeshContinuum* grid,
                  chi_math::UnknownManager* unknown_manager)
{
  unsigned int N = 1;

  if (unknown_manager != nullptr)
    N = unknown_manager->GetTotalUnknownSize();

  return grid->cells.GetNumGhosts()*N;
}

//###################################################################
/**Returns the ghost DOF indices.*/
std::vector<int> SpatialDiscretization_FV::
  GetGhostDOFIndices(chi_mesh::MeshContinuum* grid,
                     chi_math::UnknownManager* unknown_manager,
                     unsigned int unknown_id)
{
  std::vector<int> dof_ids;
  std::vector<int> ghost_cell_ids = grid->cells.GetGhostGlobalIDs();

  unsigned int N = 1;

  if (unknown_manager != nullptr)
    N = unknown_manager->GetTotalUnknownSize();

  dof_ids.reserve(GetNumGhostDOFs(grid,unknown_manager));
  for (int cell_id : ghost_cell_ids)
    for (int comp=0; comp<N; ++comp)
    {
      dof_ids.push_back(
        MapDOF(grid->cells[cell_id],unknown_manager,unknown_id,comp) );
    }

  return dof_ids;
}