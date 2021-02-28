#include "diffusion_solver.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/VolumeMesher/chi_volumemesher.h"
#include "ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h"
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include "../Boundaries/chi_diffusion_bndry_dirichlet.h"
#include "../Boundaries/chi_diffusion_bndry_reflecting.h"

#include <chi_mpi.h>
#include <chi_log.h>

extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;




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
      boundaries.push_back(new chi_diffusion::BoundaryDirichlet);
      chi_log.Log(LOG_0VERBOSE_1)
      << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";
    }

    boundaries.push_back(new chi_diffusion::BoundaryReflecting);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Reflecting boundary added (index " << boundaries.size()-1 <<  ").";

    boundaries.push_back(new chi_diffusion::BoundaryReflecting);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Reflecting boundary added (index " << boundaries.size()-1 <<  ").";
  }
  else if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
  {
    boundaries.push_back(new chi_diffusion::BoundaryDirichlet);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";

    boundaries.push_back(new chi_diffusion::BoundaryDirichlet);
    chi_log.Log(LOG_0VERBOSE_1)
      << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";
  }
  else
  {
    for (int b=0; b<std::max(1,(int)region->boundaries.size()); b++)
    {
      boundaries.push_back(new chi_diffusion::BoundaryDirichlet);
      chi_log.Log(LOG_0VERBOSE_1)
        << "Dirichlet boundary added (index " << boundaries.size()-1 <<  ").";
    }
  }

  common_items_initialized = true;
}


