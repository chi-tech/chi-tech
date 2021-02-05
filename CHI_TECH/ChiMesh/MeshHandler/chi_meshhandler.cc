#include "chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Region/chi_region.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Obtains a pointer to the last created grid.*/
chi_mesh::MeshContinuumPtr chi_mesh::MeshHandler::GetGrid(int region_index)
{
  if (region_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_mesh::MeshHandler::GetGrid. No regions added to the handler.";
    exit(EXIT_FAILURE);
  }

  if (region_index < 0)
    return region_stack.back()->GetGrid();
  else
  {
    if (region_index >= region_stack.size())
    {
      chi_log.Log(LOG_ALLERROR)
        << "Call to chi_mesh::MeshHandler::GetGrid with invalid region index "
        << region_index;
      exit(EXIT_FAILURE);
    }

    return region_stack[region_index]->GetGrid();
  }
}