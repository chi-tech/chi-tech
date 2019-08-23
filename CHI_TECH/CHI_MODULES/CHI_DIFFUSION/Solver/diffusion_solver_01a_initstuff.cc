#include "diffusion_solver.h"

#include "../../../CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../../CHI_MESH/CHI_VOLUMEMESHER/chi_volumemesher.h"
#include "../../../CHI_MESH/CHI_VOLUMEMESHER/Extruder/volmesher_extruder.h"
#include <CHI_MESH/CHI_VOLUMEMESHER/Linemesh1D/volmesher_linemesh1d.h>
#include "../Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../Boundaries/chi_diffusion_bndry_reflecting.h"

#include <chi_log.h>

extern CHI_LOG chi_log;

//###################################################################
/**Initialization of common to all solver types.*/
void chi_diffusion::Solver::InitializeCommonItems()
{
  chi_mesh::Region*  region = this->regions.back();

  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherExtruder))
  {
    for (int b=0; b<(region->boundaries.size()-2); b++)
    {
      chi_diffusion::Boundary* new_bndry =
        new chi_diffusion::BoundaryDirichlet;
      this->boundaries.push_back(new_bndry);
      chi_log.Log(LOG_0VERBOSE_1)
      << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";
    }
    chi_diffusion::Boundary* new_bndry;

    new_bndry = new chi_diffusion::BoundaryReflecting;
    this->boundaries.push_back(new_bndry);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Reflecting boundary added (index " << boundaries.size()-1 <<  ").";

    new_bndry = new chi_diffusion::BoundaryReflecting;
    this->boundaries.push_back(new_bndry);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Reflecting boundary added (index " << boundaries.size()-1 <<  ").";
  }
  else if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
  {
    chi_diffusion::Boundary* new_bndry =
      new chi_diffusion::BoundaryDirichlet;
    this->boundaries.push_back(new_bndry);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";

    new_bndry =
      new chi_diffusion::BoundaryDirichlet;
    this->boundaries.push_back(new_bndry);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";
  }
  else
  {
    for (int b=0; b<region->boundaries.size(); b++)
    {
      chi_diffusion::Boundary* new_bndry =
        new chi_diffusion::BoundaryDirichlet;
      this->boundaries.push_back(new_bndry);
      chi_log.Log(LOG_0VERBOSE_1)
        << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";
    }
  }

  common_items_initialized = true;
}