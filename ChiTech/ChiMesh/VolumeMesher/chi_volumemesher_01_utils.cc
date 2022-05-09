#include "chi_volumemesher.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/** Sets the grid member of the volume mesher.*/
void chi_mesh::VolumeMesher::SetContinuum(MeshContinuumPtr &grid)
{
  m_grid = grid;
}

//###################################################################
/** Gets the smart-pointer for the grid.*/
chi_mesh::MeshContinuumPtr& chi_mesh::VolumeMesher::GetContinuum()
{
  return m_grid;
}
