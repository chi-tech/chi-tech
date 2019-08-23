#include "surfmesher_predefined.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../CHI_REGION/chi_region.h"
#include "../../CHI_BOUNDARY/chi_boundary.h"
#include<iostream>

#include <chi_log.h>

extern CHI_LOG chi_log;

void chi_mesh::SurfaceMesherPredefined::Execute()
{
  chi_log.Log(LOG_0VERBOSE_1) << "SurfaceMesherPredefined executed";

  //================================================== Get the current handler
  chi_mesh::MeshHandler* mesh_handler = chi_mesh::GetCurrentHandler();

  //================================================== Loop over all regions
  std::vector<chi_mesh::Region*>::iterator region_iter;
  for (region_iter = mesh_handler->region_stack.begin();
       region_iter != mesh_handler->region_stack.end();
       region_iter++)
  {
    chi_mesh::Region* region = *region_iter;
    //=========================================== Check for interfaces
    //=========================================== Clear non-initial continuums
    region->volume_mesh_continua.clear();

    //=========================================== Create new continuum
    chi_mesh::MeshContinuum* remeshed_surfcont = new chi_mesh::MeshContinuum;
    region->volume_mesh_continua.push_back(remeshed_surfcont);

  }


}