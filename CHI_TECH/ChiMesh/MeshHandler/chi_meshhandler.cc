#include "chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Region/chi_region.h"

#include "chi_log.h"
extern ChiLog chi_log;

//###################################################################
/**Obtains a pointer to the last created grid.*/
chi_mesh::MeshContinuum* chi_mesh::MeshHandler::GetGrid()
{
  if (region_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_mesh::MeshHandler::GetGrid. No regions added to the handler.";
    exit(EXIT_FAILURE);
  }

  if (region_stack.back()->volume_mesh_continua.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_mesh::MeshHandler::GetGrid. Region found but no grids "
         "added to the region.";
    exit(EXIT_FAILURE);
  }

  return region_stack.back()->volume_mesh_continua.back();
}