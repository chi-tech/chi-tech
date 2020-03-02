#include "lbs_linear_boltzman_solver.h"
#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/VolumeMesher/Linemesh1D/volmesher_linemesh1d.h>
#include <ChiMesh/VolumeMesher/Extruder/volmesher_extruder.h>
#include <ChiMesh/VolumeMesher/Predefined2D/volmesher_predefined2d.h>
#include <ChiMesh/VolumeMesher/Predefined3D/volmesher_predefined3d.h>
#include <ChiMesh/SurfaceMesher/Triangle/triangle_mesher.h>
#include <ChiMesh/SurfaceMesher/Predefined/surfmesher_predefined.h>

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/** Computes the number of moments for the given mesher types*/
void LinearBoltzman::Solver::ComputeNumberOfMoments()
{
  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
  {
    int L = options.scattering_order;
    this->num_moments = L+1;
  }
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherExtruder) or
      typeid(*mesher) == typeid(chi_mesh::VolumeMesherPredefined2D) or
      typeid(*mesher) == typeid(chi_mesh::VolumeMesherPredefined3D))
  {
    int L = options.scattering_order;
    this->num_moments = L*(L+2) + 1;
  }
}

