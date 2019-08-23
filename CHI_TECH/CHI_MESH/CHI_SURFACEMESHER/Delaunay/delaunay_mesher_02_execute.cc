#include "delaunay_mesher.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../CHI_REGION/chi_region.h"
#include "../../../CHI_TIMER/chi_timer.h"
#include "../../CHI_BOUNDARY/chi_boundary.h"
#include "../../CHI_MESHCONTINUUM/chi_meshcontinuum.h"


/**Executes a Delaunay surface remeshing operation on all
 * regions.
 *
 * ## Concepts
 *  - Create a Box of contexts for each region (RegionDelaunyContexts)
 *      - Add a DelaunayMeshContext for each boundary
 *  - Hash out interfaces between regions (Currently on hold)
 *  - Mesh first RegionDelaunyContexts
 *      - For each surface mesh
 *          - Split into into patches (Co-planar and contigous)
 *          - For each patch
 *              - Get edge loops and essential vertices
 *              - If not first patch
 *                  - Add vertices from in-region interfaces
 *              - Project to 2D
 *              - Perform Lexicographical meshing
 * */
void chi_mesh::SurfaceMesherDelaunay::Execute()
{
  std::cout << "SurfaceMesherDelaunay executed";
  std::cout << std::endl;

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
    chi_mesh::Region* region = *region_iter;

    //=========================================== Create RegionContext
    RegionDelaunyContexts* region_context = new RegionDelaunyContexts;
    region_contexts.push_back(region_context);

    //=========================================== Check for interfaces
    //=========================================== Clear non-initial continuums
    //region->mesh_continua.clear();

    //=========================================== Create new continuum
    //chi_mesh::MeshContinuum* remeshed_surfcont = new chi_mesh::MeshContinuum;
    //region->mesh_continua.push_back(remeshed_surfcont);

    //=========================================== Create Delaunay context for each boundary
    std::vector<chi_mesh::Boundary*>::iterator bndry;
    for (bndry = region->boundaries.begin();
         bndry != region->boundaries.end();
         bndry++)
    {
      DelaunayMeshContext* new_del_context
        = new DelaunayMeshContext;
      new_del_context->context_boundary = (*bndry);
      if ((*bndry)->initial_mesh_continuum.line_mesh!= nullptr)
      {
        new_del_context->context_line_mesh = (*bndry)->initial_mesh_continuum.line_mesh;
      }
      if ((*bndry)->initial_mesh_continuum.surface_mesh!= nullptr)
      {
        new_del_context->context_surface_mesh = (*bndry)->initial_mesh_continuum.surface_mesh;
      }

      region_context->push_back(new_del_context);
    }
  }

  //================================================== Loop over regioncontexts
  //Mesh each region
  std::vector<RegionDelaunyContexts*>::iterator region_context;
  for (region_context = region_contexts.begin();
       region_context != region_contexts.end();
       region_context++)
  {
    printf("Meshing region %d\n",std::distance(region_contexts.begin(),region_context));
    MeshRegion(*(*region_context));
  }

  //================================================== Copy surface meshes
  //============================================= Loop over regions
  for (unsigned r=0; r<mesh_handler->region_stack.size(); r++)
  {
    chi_mesh::Region*          cur_region = mesh_handler->region_stack[r];
    RegionDelaunyContexts* cur_region_context = region_contexts[r];

    //====================================== Loop over boundaries
    for (unsigned b=0; b<cur_region->boundaries.size(); b++)
    {
      if (cur_region->boundaries[b]->initial_mesh_continuum.surface_mesh != nullptr)
      {
        //============================= Clear the mesh continua
        cur_region->boundaries[b]->mesh_continua.clear();

        //============================= Create and assign new mesh continua
        MeshContinuum* new_mesh_cont = new MeshContinuum;
        new_mesh_cont->surface_mesh = cur_region_context->at(b)->remeshed_surface_mesh;
        printf("Surface mesh copied with %d triangles and %d vertices\n", new_mesh_cont->surface_mesh->faces.size(),
               new_mesh_cont->surface_mesh->vertices.size());
        cur_region->boundaries[b]->mesh_continua.push_back(new_mesh_cont);
      }
    }
  }



  //======================================================= Getting end-time
  printf("SurfaceMesherDelaunay execution completed. %.4f s\n",tmesh_timer.GetTime()/1000.0);


}