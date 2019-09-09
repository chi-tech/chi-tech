#include"../../../ChiLua/chi_lua.h"
#include <iostream>
#include "../chi_region.h"
#include "../../CHI_SURFACEMESH/chi_surfacemesh.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../CHI_BOUNDARY/chi_boundary.h"



//#############################################################################
/** Obtains a handle to the surface mesh associated with a boundary

\param RegionHandle int Handle to the region for which boundary is to be added.
\param BoundaryNumber int Index of the boundary for which the mesh is to be obtained.
\param ContinuumNumber int Optional. If supplied then the surface mesh for the given con
 tinuum will be used

\return SurfaceMeshHandle int. Handle to the surface mesh extracted.

\ingroup LuaRegion
\author Jan*/
int chiRegionGetBoundarySurfaceMesh(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int num_args = lua_gettop(L);
  int region_index = lua_tonumber(L,1);
  int boundary_index = lua_tonumber(L,2);

  chi_mesh::Region* cur_region;
  chi_mesh::Boundary* cur_boundary;
  try{
    cur_region = cur_hndlr->region_stack.at(region_index);
  }
  catch(const std::invalid_argument& ia)
  {
    std::cerr << "ERROR: Invalid index to region.\n";
    exit(EXIT_FAILURE);
  }

  try{
    cur_boundary = cur_region->boundaries.at(boundary_index);
  }
  catch(const std::invalid_argument& ia)
  {
    std::cerr << "ERROR: Invalid index to boundary.\n";
    exit(EXIT_FAILURE);
  }

  if (num_args == 2)
  {
    if (cur_boundary->initial_mesh_continuum.surface_mesh != nullptr)
    {
      chi_mesh::SurfaceMesh* cur_surfmesh = cur_boundary->initial_mesh_continuum.surface_mesh;
      cur_hndlr->surface_mesh_stack.push_back(cur_surfmesh);

      lua_pushnumber(L,cur_hndlr->surface_mesh_stack.size()-1);
      printf("Hello\n");
      return 1;
    } else
    {
      std::cerr << "ERROR: No surface mesh associated with boundary.\n";
      exit(EXIT_FAILURE);
    }

  }
  else if (num_args == 3)
  {
    int continuum_index = lua_tonumber(L,3);
    chi_mesh::MeshContinuum* cur_cont;
    try{
      cur_cont = cur_boundary->mesh_continua.at(continuum_index);
    }
    catch(const std::invalid_argument& ia)
    {
      std::cerr << "ERROR: Invalid index to mesh continuum.\n";
      exit(EXIT_FAILURE);
    }
    if (cur_cont->surface_mesh != nullptr)
    {
      chi_mesh::SurfaceMesh* cur_surfmesh = cur_cont->surface_mesh;
      cur_hndlr->surface_mesh_stack.push_back(cur_surfmesh);

      lua_pushnumber(L,cur_hndlr->surface_mesh_stack.size()-1);
      return 1;
    } else
    {
      std::cerr << "ERROR: No surface mesh associated with boundary.\n";
      exit(EXIT_FAILURE);
    }
  }



  return 0;
}