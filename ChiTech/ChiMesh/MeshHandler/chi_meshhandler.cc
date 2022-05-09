#include "chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Obtains a pointer to the last created grid.*/
chi_mesh::MeshContinuumPtr& chi_mesh::MeshHandler::GetGrid() const
{
  if (volume_mesher == nullptr)
    throw std::logic_error("chi_mesh::MeshHandler::GetGrid: Volume mesher "
                           "undefined. Requires a volume mesher to be "
                           "created.");

  return volume_mesher->GetContinuum();
}