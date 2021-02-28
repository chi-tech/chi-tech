#include "surfmesher_passthrough.h"
#include "../../MeshHandler/chi_meshhandler.h"
#include "../../Region/chi_region.h"
#include "../../Boundary/chi_boundary.h"
#include<iostream>

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/**Executes the pass-through surface mesher.*/
void chi_mesh::SurfaceMesherPassthrough::Execute()
{
  chi_log.Log(LOG_0VERBOSE_1) << "SurfaceMesherPassthrough executed";

  //================================================== Get the current handler
  chi_mesh::MeshHandler* mesh_handler = chi_mesh::GetCurrentHandler();

  //================================================== Check empty region list
  if (mesh_handler->region_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "SurfaceMesherPassthrough: No region added.";
    exit(EXIT_FAILURE);
  }

  //================================================== Loop over all regions
  for (auto region : mesh_handler->region_stack)
  {
  }

}