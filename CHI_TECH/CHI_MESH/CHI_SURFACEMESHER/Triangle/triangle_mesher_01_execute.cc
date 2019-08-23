#include "triangle_mesher.h"
#include "../../../CHI_TIMER/chi_timer.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../CHI_REGION/chi_region.h"
#include "../../CHI_BOUNDARY/chi_boundary.h"
#include "../../../CHI_MPI/chi_mpi.h"
#include <chi_log.h>

extern CHI_MPI chi_mpi;
extern CHI_LOG chi_log;

/**Executes the triangle mesher*/
void chi_mesh::SurfaceMesherTriangle::Execute()
{
  chi_log.Log(LOG_0) << "SurfaceMesherTriangle executed";

  //================================================== Start timer
  CHI_TIMER tmesh_timer;
  tmesh_timer.Reset();

  //================================================== Get the current handler
  chi_mesh::MeshHandler* mesh_handler = chi_mesh::GetCurrentHandler();

  //================================================== Loop over all regions
  std::vector<chi_mesh::Region*>::iterator region_iter;
  for (region_iter = mesh_handler->region_stack.begin();
       region_iter != mesh_handler->region_stack.end();
       region_iter++)
  {
    chi_mesh::Region *region = *region_iter;

    //=========================================== Loop over boundary
    std::vector<chi_mesh::Boundary*>::iterator bndry;
    for (bndry = region->boundaries.begin();
         bndry != region->boundaries.end();
         bndry++)
    {
      if ((*bndry)->initial_mesh_continuum.surface_mesh != nullptr)
      {

        if (this->get_auto_min)
        {CalculateMinimumArea((*bndry));}

        AddCutLines((*bndry));

        //(*bndry)->initial_mesh_continuum.surface_mesh->ExportToOBJFile("Test.obj");

        chi_log.Log(LOG_ALLVERBOSE_2) << "Triangulating...";
        MeshBoundary((*bndry));
      }
    }
  }



  //======================================================= Getting end-time
  chi_log.Log(LOG_0) << "SurfaceMesherTriangle execution "
                                         "completed. "
                                      << tmesh_timer.GetTime()/1000.0
                                      << " s";
  MPI_Barrier(MPI_COMM_WORLD);

}