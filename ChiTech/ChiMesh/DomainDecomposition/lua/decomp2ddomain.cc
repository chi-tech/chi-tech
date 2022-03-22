#include "ChiLua/chi_lua.h"
#include<iostream>
#include "ChiMesh/MeshHandler/chi_meshhandler.h"


/** \defgroup LuaDomainDecomposition Domain decomposition
 * \ingroup LuaMesh
*/

//#############################################################################
/** Decomposes a surface mesh into block px py elements.
 * \image html "InProgressImage.png" width=200px
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