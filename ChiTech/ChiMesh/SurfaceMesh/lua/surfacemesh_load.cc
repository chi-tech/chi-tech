#include "ChiLua/chi_lua.h"

#include <iostream>
#include "ChiMesh/SurfaceMesh/chi_surfacemesh.h"
#include "../../MeshHandler/chi_meshhandler.h"

#include "chi_runtime.h"

#include "chi_log.h"
extern ChiLog& chi_log;


//############################################################################# Create Window
/** Loads mesh data from a wavefront object.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.
\param FileName char* Path to the file to be imported.
\param polyflag bool (Optional)Flag indicating whether triangles
 are to be read as polygons. [Default: true)

\return success bool Return true if file was successfully loaded and false
 otherwise.
\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshImportFromOBJFile(lua_State *L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Get arguments
  int num_args = lua_gettop(L);
  int handle = lua_tonumber(L,1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  bool as_poly = true;
  if (num_args==3)
  {
    as_poly = lua_toboolean(L,3);
  }

  auto& surface_mesh = chi::GetStackItem<chi_mesh::SurfaceMesh>(
    chi::surface_mesh_stack, handle, __FUNCTION__);

  std::stringstream outtext;
  outtext << "chiSurfaceMeshImportFromOBJFile: "
             "Loading Wavefront .obj file: ";
  outtext << temp << std::endl;
  chi_log.Log(LOG_ALLVERBOSE_2) << outtext.str();
  surface_mesh.ImportFromOBJFile(temp, as_poly);

  return 1;
}

//############################################################################# Create Window
/** Loads mesh data from a wavefront object.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.
\param FileName char* Path to the file to be imported.
\param polyflag bool (Optional)Flag indicating whether triangles
 are to be read as polygons. [Default: true)

\return success bool Return true if file was successfully loaded and false
 otherwise.
\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshImportFromTriangleFiles(lua_State *L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Get arguments
  int num_args = lua_gettop(L);
  int handle = lua_tonumber(L,1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  bool as_poly = true;
  if (num_args==3)
  {
    as_poly = lua_toboolean(L,3);
  }

  auto& surface_mesh = chi::GetStackItem<chi_mesh::SurfaceMesh>(
    chi::surface_mesh_stack, handle, __FUNCTION__);

  surface_mesh.ImportFromTriangleFiles(temp,as_poly);

  return 1;
}

int chiSurfaceMeshImportFromMshFiles(lua_State *L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Get arguments
  int num_args = lua_gettop(L);
  int handle = lua_tonumber(L,1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  bool as_poly = true;
  if (num_args==3)
  {
    as_poly = lua_toboolean(L,3);
  }

  auto& surface_mesh = chi::GetStackItem<chi_mesh::SurfaceMesh>(
    chi::surface_mesh_stack, handle, __FUNCTION__);

  std::stringstream outtext;
  outtext << "chiSurfaceMeshImportFromMshFiles: "
             "Loading a gmsh ascii file: ";
  outtext << temp << std::endl;
  chi_log.Log(LOG_ALLVERBOSE_2) << outtext.str();
  surface_mesh.ImportFromMshFiles(temp, as_poly);

  return 1;
}
