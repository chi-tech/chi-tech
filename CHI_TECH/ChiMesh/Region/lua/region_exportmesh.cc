#include"../../../ChiLua/chi_lua.h"
#include <iostream>
#include "../chi_region.h"
#include "../../SurfaceMesh/chi_surfacemesh.h"
#include "../../MeshHandler/chi_meshhandler.h"
#include "../../Boundary/chi_boundary.h"

#include <chi_log.h>
extern ChiLog chi_log;

//#############################################################################
/** Exports the mesh to python.

\param RegionHandle int Handle to the region for which boundary is to be added.
\param FileName char Name of the file to be used.
\param ExportTemplate bool Default: False. Flag indicating whether to export
                     the extruder's surface mesh template.

\ingroup LuaRegion
\author Jan*/
int chiRegionExportMeshToPython(lua_State *L)
{
//  //============================================= Check arguments
//  int num_args = lua_gettop(L);
//  if (!((num_args == 2) || (num_args == 3)))
//  {
//    chi_log.Log(LOG_0ERROR) << "Incorrect amount of arguments used in "
//                               "chiRegionExportMeshToPython";
//    exit(EXIT_FAILURE);
//  }
//
//  int region_index = lua_tonumber(L,1);
//  const char* file_name = lua_tostring(L,2);
//  bool export_template = false;
//
//  if (num_args == 3) export_template = lua_toboolean(L,3);
//
//  //============================================= Get current handler
//  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();
//
//  //============================================= Attempt to obtain region
//  chi_mesh::Region* cur_region;
//  try{
//    cur_region = cur_hndlr->region_stack.at(region_index);
//  }
//  catch(const std::invalid_argument& ia)
//  {
//    chi_log.Log(LOG_0ERROR) << "ERROR: Invalid index to region in "
//                               "chiRegionExportMeshToPython.";
//    exit(EXIT_FAILURE);
//  }
//
//  //============================================= Get back continuum
//  if (cur_region->volume_mesh_continua.size()>0)
//  {
//    int num_cont = cur_region->volume_mesh_continua.size();
//
//    chi_mesh::MeshContinuum* vol_cont;
//    if ((export_template) && (num_cont >= 2))
//      vol_cont = cur_region->volume_mesh_continua[num_cont-2];
//    else
//      vol_cont= cur_region->volume_mesh_continua.back();
//
//    vol_cont->ExportCellsToPython((char*)file_name);
//  }
//  else
//  {
//    chi_log.Log(LOG_ALLWARNING) << "No volume continuum to export in "
//                                   "call to chiRegionExportMeshToPython.";
//  }
//
  chi_log.Log(LOG_0WARNING)
    << "chiRegionExportMeshToPython is deprecated. "
       "Use chiRegionExportMeshToVTK instead.";

  return 0;
}



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
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Attempt to obtain region
  chi_mesh::Region* cur_region;
  try{
    cur_region = cur_hndlr->region_stack.at(region_index);
  }
  catch(const std::invalid_argument& ia)
  {
    chi_log.Log(LOG_0ERROR) << "ERROR: Invalid index to region in "
                               "chiRegionExportMeshToObj.";
    exit(EXIT_FAILURE);
  }

  auto vol_cont= cur_region->GetGrid();
  vol_cont->ExportCellsToObj((char*)file_name,per_material);


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
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Attempt to obtain region
  chi_mesh::Region* cur_region;
  try{
    cur_region = cur_hndlr->region_stack.at(region_index);
  }
  catch(const std::invalid_argument& ia)
  {
    chi_log.Log(LOG_0ERROR) << "ERROR: Invalid index to region in "
                               "chiRegionExportMeshToObj.";
    exit(EXIT_FAILURE);
  }

  auto vol_cont = cur_region->GetGrid();

  vol_cont->ExportCellsToVTK(base_name);

  return 0;
}