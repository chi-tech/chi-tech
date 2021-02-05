#include <ChiLua/chi_lua.h>
#include "../surfacemesher.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>
#include <ChiMesh/Region/chi_region.h>
#include <ChiMesh/Boundary/chi_boundary.h>
#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/**Exports the first available surface mesh to a wavefront .obj file.
 * This would be the surface mesh associated with the last mesh operation.
\ingroup LuaSurfaceMesher
 * */
int chiSurfaceMesherExportToObj(lua_State* L)
{
  int num_args = lua_gettop(L);

  //================================================== Extract file name
  const char* file_name = lua_tostring(L,1);

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
      if ((*bndry)->initial_mesh_continuum->surface_mesh != nullptr)
      {
        chi_mesh::MeshContinuum* twod_grid = (*bndry)->mesh_continua.back();

        twod_grid->surface_mesh->ExportToOBJFile(file_name);
      }
    }
  }

  return 0;
}
