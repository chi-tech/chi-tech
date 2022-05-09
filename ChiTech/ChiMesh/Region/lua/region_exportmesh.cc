#include "ChiLua/chi_lua.h"
#include "ChiMesh/SurfaceMesh/chi_surfacemesh.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//#############################################################################
/** Exports the mesh to obj format.

\param RegionHandle int Handle to the region for which boundary is to be added.
\param FileName char Name of the file to be used.
\param ExportByMaterial bool Default: False. Flag indicating whether to export
                     the extruder's surface mesh by material.

\ingroup LuaRegion
\author Jan*/
int chiRegionExportMeshToObj(lua_State *L)
{
  //============================================= Check arguments
  int num_args = lua_gettop(L);
  if (!((num_args == 2) || (num_args == 3)))
  {
    chi_log.Log(LOG_0ERROR) << "Incorrect amount of arguments used in "
                               "chiRegionExportMeshToObj";
    exit(EXIT_FAILURE);
  }

  int region_index = lua_tonumber(L,1);
  const char* file_name = lua_tostring(L,2);
  bool per_material = false;

  if (num_args == 3) per_material = lua_toboolean(L,3);

  //============================================= Get current handler
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  auto& grid = cur_hndlr.GetGrid();
  grid->ExportCellsToObj((char*)file_name, per_material);


  return 0;
}


//#############################################################################
/** Exports the mesh to vtu format.

\param RegionHandle int Handle to the region for which boundary is to be added.
\param FileName char Name of the file to be used.

\ingroup LuaRegion
\author Jan*/
int chiRegionExportMeshToVTK(lua_State *L)
{
  //============================================= Check arguments
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("chiRegionExportMeshToVTK",2,num_args);

  int region_index = lua_tonumber(L,1);
  const char* base_name = lua_tostring(L,2);


  //============================================= Get current handler
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  auto& grid = cur_hndlr.GetGrid();

  grid->ExportCellsToVTK(base_name);

  return 0;
}