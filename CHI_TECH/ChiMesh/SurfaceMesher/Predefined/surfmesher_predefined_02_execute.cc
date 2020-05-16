#include "surfmesher_predefined.h"
#include "../../MeshHandler/chi_meshhandler.h"
#include "../../Region/chi_region.h"
#include "../../Boundary/chi_boundary.h"
#include<iostream>

#include <chi_log.h>

extern ChiLog& chi_log;

void chi_mesh::SurfaceMesherPredefined::Execute()
{
  chi_log.Log(LOG_0VERBOSE_1) << "SurfaceMesherPredefined executed";

  //================================================== Get the current handler
  chi_mesh::MeshHandler* mesh_handler = chi_mesh::GetCurrentHandler();

  //================================================== Check empty region list
  if (mesh_handler->region_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "SurfaceMesherPredefined: No region added.";
    exit(EXIT_FAILURE);
  }

  //================================================== Loop over all regions
//  std::vector<chi_mesh::Region*>::iterator region_iter;
//  for (region_iter = mesh_handler->region_stack.begin();
//       region_iter != mesh_handler->region_stack.end();
//       region_iter++)
  for (auto region : mesh_handler->region_stack)
  {
//    chi_mesh::Region* region = *region_iter;
    //=========================================== Check for interfaces
    //=========================================== Clear non-initial continuums
//    region->volume_mesh_continua.clear();

    //=========================================== Create new continuum
//    chi_mesh::MeshContinuum* remeshed_surfcont = new chi_mesh::MeshContinuum;
//    region->volume_mesh_continua.push_back(remeshed_surfcont);

  }


}