#include "ChiLua/chi_lua.h"

#include "../chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Exports the mesh to a wavefront.obj format.
\param FileName char Name of the file to be used.
\param ExportByMaterial bool Default: False. Flag indicating whether to export
                     the extruder's surface mesh by material.
\ingroup LuaMeshHandler
*/
int chiMeshHandlerExportMeshToObj(lua_State* L)
{
  //============================================= Check arguments
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(fname, 1, num_args);

  const std::string file_name = lua_tostring(L,1);

  bool per_material = false;
  if (num_args == 2) per_material = lua_toboolean(L,2);

  //============================================= Get current handler
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  auto& grid = cur_hndlr.GetGrid();
  grid->ExportCellsToObj(file_name.c_str(), per_material);

  return 0;
}


//###################################################################
/**Exports the mesh to a wavefront.obj format.
\param FileName char Name of the file to be used.
\ingroup LuaMeshHandler
*/
int chiMeshHandlerExportMeshToVTK(lua_State* L)
{
  //============================================= Check arguments
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  const std::string file_name = lua_tostring(L,1);

  //============================================= Get current handler
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  auto& grid = cur_hndlr.GetGrid();
  grid->ExportCellsToVTK(file_name.c_str());

  return 0;
}