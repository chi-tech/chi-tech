#include "chi_nptransport.h"
#include <CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h>
#include <CHI_MESH/CHI_VOLUMEMESHER/Linemesh1D/volmesher_linemesh1d.h>
#include <CHI_MESH/CHI_VOLUMEMESHER/Extruder/volmesher_extruder.h>
#include <CHI_MESH/CHI_VOLUMEMESHER/Predefined2D/volmesher_predefined2d.h>
#include <CHI_MESH/CHI_SURFACEMESHER/Triangle/triangle_mesher.h>
#include <CHI_MESH/CHI_SURFACEMESHER/Predefined/surfmesher_predefined.h>

#include <chi_log.h>

extern CHI_LOG chi_log;

//###################################################################
/** Computes the number of moments for the given mesher types*/
void CHI_NPTRANSPORT::ComputeNumberOfMoments()
{
  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;

  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherLinemesh1D))
  {
    int L = options.scattering_order;
    this->num_moments = L+1;
  }
  if (typeid(*mesher) == typeid(chi_mesh::VolumeMesherExtruder) or
      typeid(*mesher) == typeid(chi_mesh::VolumeMesherPredefined2D))
  {
    int L = options.scattering_order;
    this->num_moments = L*(L+2) + 1;
  }
}

//###################################################################
/**Sets partitioning values.*/
void CHI_NPTRANSPORT::SetPartitioning()
{
  chi_mesh::Region*              aregion = this->regions.back();
  chi_mesh::MeshContinuum* vol_continuum;

  chi_mesh::MeshHandler*    mesh_handler = chi_mesh::GetCurrentHandler();
  chi_mesh::VolumeMesher*         mesher = mesh_handler->volume_mesher;
  chi_mesh::SurfaceMesher*   surf_mesher = mesh_handler->surface_mesher;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Serial
  if (options.partition_method == PARTITION_METHOD_SERIAL)
  {

    //mesher->Execute();
    //printf("Hello Serial %d\n", aregion->volume_mesh_continua.size());
    vol_continuum = aregion->volume_mesh_continua.back();
    for (int i=0; i<vol_continuum->cells.size(); i++)
    {
      local_cell_indices.push_back(i);
    }
   // printf("Bye Serial\n");
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% From surface mesh
  else if (options.partition_method == PARTITION_METHOD_FROM_SURFACE)
  {
    //printf("Hello partition method from surface\n");
    if ((typeid(*surf_mesher) != typeid(chi_mesh::SurfaceMesherTriangle)) &&
        (typeid(*surf_mesher) != typeid(chi_mesh::SurfaceMesherPredefined)))
    {
      chi_log.Log(LOG_ALLERROR) <<
                     "PARTITION_METHOD_FROM_SURFACE can only be used"
                     "with SURFACEMESHER_TRIANGLE or SURFACEMESHER_PREDEFINED.";
      exit(EXIT_FAILURE);
    }

    if ((typeid(*mesher) != typeid(chi_mesh::VolumeMesherExtruder)) &&
        (typeid(*mesher) != typeid(chi_mesh::VolumeMesherPredefined2D)))
    {
      chi_log.Log(LOG_ALLERROR) <<
                     "PARTITION_METHOD_FROM_SURFACE can only be used"
                     "with SURFACEMESHER_TRIANGLE or "
                     "VOLUMEMESHER_PREDEFINED2D.";
      exit(EXIT_FAILURE);
    }

    //mesher->Execute();
    vol_continuum = aregion->volume_mesh_continua.back();
    //================================================ Copy cell indices
    //                                                 from mesher
//    for (int i=0; i<mesher->local_cell_glob_indices.size(); i++)
//    {
//      local_cell_indices.push_back(mesher->local_cell_glob_indices[i]);
//    }
  }
  else
  {
    fprintf(stderr,"Unknown partition method specified.\n");
    exit(EXIT_FAILURE);
  }
}