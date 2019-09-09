#include"../../../ChiLua/chi_lua.h"
#include<iostream>
#include "../chi_domdecomp.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"
#include "../../CHI_SURFACEMESHER/Triangle/triangle_mesher.h"

#include <typeinfo>

/** \defgroup LuaDomainDecomposition Domain decomposition
 * \ingroup LuaMesh
*/

//#############################################################################
/** Decomposes a region domain and stores a collection of lists that indicates
 * which cells go to which process.
 *
\param Px int Number of divisions in x.
\param Py int Number of divisions in y.
\param RegionHandle int Handle to the region that is to be decomposed.

\return Handle int Handle to the created collecion.
\ingroup LuaDomainDecomposition
\author Jan*/
int chiDomDecompose2D(lua_State *L)
{
  //================================================== Get current handler
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //================================================== Extract arguments
  int px           = lua_tonumber(L,1);
  int py           = lua_tonumber(L,2);
  int region_index = lua_tonumber(L,3);

  chi_mesh::Region* cur_region;
  try{
    cur_region = cur_hndlr->region_stack.at(region_index);
  }
  catch(const std::invalid_argument& ia)
  {
    std::cerr << "ERROR: Invalid index to region.\n";
    exit(EXIT_FAILURE);
  }

  fprintf(stdout,"Decomposing region %d\n",region_index);

  //================================================== Get surface mesher
  if (typeid(*cur_hndlr->surface_mesher) ==
      typeid(chi_mesh::SurfaceMesherTriangle))
  {
    chi_mesh::Decompose2DDomain(px,py,
      (chi_mesh::SurfaceMesherTriangle*)cur_hndlr->surface_mesher,
      cur_region);
  }
  else
  {
    fprintf(stderr,"Invalid mesher used to perform domain decompisition\n");
    exit(EXIT_FAILURE);
  }






  fprintf(stdout,"Domain decomposition completed, collection index %d\n",region_index);

  return 1;
}

//#############################################################################
/** Decomposes a surface mesh into block px py elements.
 *
\param Surface mesh handler
\param Px int Number of divisions in x.
\param Py int Number of divisions in y.

\ingroup LuaDomainDecomposition
\author Jan*/
int chiDecomposeSurfaceMeshPxPy(lua_State *L)
{
  int num_args = lua_gettop(L);

  if (num_args != 3)
    LuaPostArgAmountError("chiDecomposeSurfaceMeshPxPy",3,num_args);

  //================================================== Get current handler
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //================================================== Extract arguments
  int surface_hndl = lua_tonumber(L,1);
  int px           = lua_tonumber(L,2);
  int py           = lua_tonumber(L,3);


  chi_mesh::SurfaceMesh* surf_mesh;
  try{
    surf_mesh = cur_hndlr->surface_mesh_stack.at(surface_hndl);
  }
  catch(const std::invalid_argument& ia)
  {
    std::cerr << "ERROR: Invalid index to surface mesh.\n";
    exit(EXIT_FAILURE);
  }

  chi_mesh::DecomposeSurfaceMeshPxPy(surf_mesh,px,py);

  return 0;
}