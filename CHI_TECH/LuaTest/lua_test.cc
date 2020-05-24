#include <ChiLua/chi_lua.h>

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//###################################################################
/**This is a lua test function.
\param argument1 Any Argument of any type.
\ingroup LuaGeneralUtilities
 */
int chiLuaTest(lua_State* L)
{
  auto mesh_handler = chi_mesh::GetCurrentHandler();
  auto grid = mesh_handler->GetGrid();

  for (auto& cell : grid->local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (face.neighbor >= 0)
        face.GetNeighborAssociatedFace(grid);
    }
  }

  return 0;
}