#include"../../../ChiLua/chi_lua.h"
#include <iostream>
#include "../chi_region.h"
#include "../../CHI_SURFACEMESH/chi_surfacemesh.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../CHI_BOUNDARY/chi_boundary.h"



//#############################################################################
/** Adds a surface mesh boundary to the region

\param RegionHandle int Handle to the region for which boundary is to be added.
\param SurfaceHandle int Handle to the surface mesh.

\return BoundaryNumber int. Number of the boundary added to the region.

\ingroup LuaRegion
\author Jan*/
int chiRegionAddSurfaceBoundary(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int region_index = lua_tonumber(L,1);
  int surface_index = lua_tonumber(L,2);

  chi_mesh::Region* cur_region;
  chi_mesh::SurfaceMesh* cur_surfacemesh;
  try{
    cur_region = cur_hndlr->region_stack.at(region_index);
  }
  catch(std::out_of_range o)
  {
    std::cerr << "ERROR: Invalid index to region.\n";
    exit(EXIT_FAILURE);
  }

  try{
    cur_surfacemesh = cur_hndlr->surface_mesh_stack.at(surface_index);
  }
  catch(std::out_of_range o)
  {
    std::cerr << "ERROR: Invalid index to surface mesh.\n";
    exit(EXIT_FAILURE);
  }

  chi_mesh::Boundary* new_boundary = new chi_mesh::Boundary;
  new_boundary->initial_mesh_continuum.surface_mesh = cur_surfacemesh;

  cur_region->boundaries.push_back(new_boundary);

  lua_pushnumber(L,cur_region->boundaries.size()-1);

  return 1;
}

//#############################################################################
/** Adds a line mesh boundary to the region

\param RegionHandle int Handle to the region for which boundary is to be added.
\param LineMeshHandle int Handle to the line mesh.

\return BoundaryNumber int. Number of the boundary added to the region.

\ingroup LuaRegion
\author Jan*/
int chiRegionAddLineBoundary(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int region_index = lua_tonumber(L,1);
  int line_index = lua_tonumber(L,2);

  chi_mesh::Region* cur_region;
  chi_mesh::LineMesh* cur_linemesh;
  try{
    cur_region = cur_hndlr->region_stack.at(region_index);
  }
  catch(std::out_of_range o)
  {
    std::cerr << "ERROR: Invalid index to region.\n";
    exit(EXIT_FAILURE);
  }

  try{
    cur_linemesh = cur_hndlr->linemesh_stack.at(line_index);
  }
  catch(std::out_of_range o)
  {
    std::cerr << "ERROR: Invalid index to line mesh.\n";
    exit(EXIT_FAILURE);
  }

  chi_mesh::Boundary* new_boundary = new chi_mesh::Boundary;
  new_boundary->initial_mesh_continuum.line_mesh = cur_linemesh;

  cur_region->boundaries.push_back(new_boundary);

  lua_pushnumber(L,cur_region->boundaries.size()-1);

  return 1;
}