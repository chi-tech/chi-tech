#include"../../../ChiLua/chi_lua.h"

#include <iostream>
#include <algorithm>
#include "../chi_surfacemesh.h"
#include "../../MeshHandler/chi_meshhandler.h"

#include <chi_log.h>

extern ChiLog& chi_log;


//#############################################################################
/** Builds sweep ordering for a number of angles and checks whether any
 * cyclic dependencies are encountered.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.
\param NumAngles int Number of azimuthal angles to use for checking cycles.

\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshCheckCycles(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("chiSurfaceMeshCheckCycles",2,num_args);

  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  int surf_handle = lua_tonumber(L,1);
  int num_angles  = lua_tonumber(L,2);

  try{
    chi_mesh::SurfaceMesh* curItem =
      cur_hndlr.surface_mesh_stack.at(surf_handle);

    curItem->CheckCyclicDependencies(num_angles);
  }

  catch(const std::out_of_range& o){
    std::cerr << "ERROR: Invalid index to surface mesh.\n";
    exit(EXIT_FAILURE);
  }
  return 0;
}

//#############################################################################
/** Computes load balancing parameters for given predictive x and y cuts
 * without actually performing cuts.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.
\param Xcuts table Array of x-values associated with the xcuts.
\param Ycuts table Array of y-values associated with the ycuts.

\ingroup LuaSurfaceMesh
\author Jan*/
int chiComputeLoadBalancing(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("chiComputeLoadBalancing",3,num_args);

  //======================================== Get reference surface mesh
  int surf_handle = lua_tonumber(L,1);
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();
  chi_mesh::SurfaceMesh* cur_surf;
  try{
    cur_surf = cur_hndlr.surface_mesh_stack.at(surf_handle);
  }
  catch(const std::out_of_range& o){
    std::cerr << "chiComputeLoadBalancing: Invalid index to surface mesh.\n";
    exit(EXIT_FAILURE);
  }

  //======================================== Extract x-cuts
  if (!lua_istable(L,2))
  {
    chi_log.Log(LOG_ALLERROR)
      << "In call to chiComputeLoadBalancing: "
      << " expected table for argument 2. Incompatible value supplied.";
    exit(EXIT_FAILURE);
  }

  int x_table_len = lua_rawlen(L,2);

  std::vector<double> x_cuts(x_table_len,0.0);
  for (int g=0; g<x_table_len; g++)
  {
    lua_pushnumber(L,g+1);
    lua_gettable(L,2);
    x_cuts[g] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  //======================================== Extract y-cuts
  if (!lua_istable(L,3))
  {
    chi_log.Log(LOG_ALLERROR)
      << "In call to chiComputeLoadBalancing: "
      << " expected table for argument 3. Incompatible value supplied.";
    exit(EXIT_FAILURE);
  }

  int y_table_len = lua_rawlen(L,3);

  std::vector<double> y_cuts(y_table_len,0.0);
  for (int g=0; g<y_table_len; g++)
  {
    lua_pushnumber(L,g+1);
    lua_gettable(L,3);
    y_cuts[g] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  //======================================== Call compute balance
  std::stable_sort(x_cuts.begin(),x_cuts.end());
  std::stable_sort(y_cuts.begin(),y_cuts.end());
  cur_surf->ComputeLoadBalancing(x_cuts,y_cuts);


  return 0;
}
