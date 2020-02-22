#include "chi_region.h"

#include "chi_log.h"

extern ChiLog chi_log;

//###################################################################
/** Obtains the latest created grid from the region.*/
chi_mesh::MeshContinuum* chi_mesh::Region::GetGrid()
{
  if (this->volume_mesh_continua.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "Region: Grid retrieval failed. Verify that volume mesher"
         " has executed for this region.";
    exit(EXIT_FAILURE);
  } else
    return volume_mesh_continua.back();

}