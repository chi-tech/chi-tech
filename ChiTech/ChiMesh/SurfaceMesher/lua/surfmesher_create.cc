#include "../../../ChiLua/chi_lua.h"
#include <iostream>
#include "../Predefined/surfmesher_predefined.h"

#include "../../MeshHandler/chi_meshhandler.h"

#include <chi_log.h>

extern ChiLog& chi_log;

//#############################################################################
/** Creates a surface preprocessor.
 *
\param SurfaceMesherType int Surface Remesher type. See SurfaceMesherType.

## _

###SurfaceMesherType:\n
SURFACEMESHER_PREDEFINED\n
 Makes no modification to the region surfaces.\n\n

\code
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
\endcode

## _

\ingroup LuaSurfaceMesher
\author Jan*/
int chiSurfaceMesherCreate(lua_State *L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Get argument
  LuaCheckNilValue("chiSurfaceMesherCreate",L,1);
  int type = lua_tonumber(L,1);

  //============================================= Create the surface mesher
  std::shared_ptr<chi_mesh::SurfaceMesher> new_mesher = nullptr;
  if (type==(int)chi_mesh::SurfaceMesherType::Predefined)
  {
    new_mesher = std::make_shared<chi_mesh::SurfaceMesherPredefined>();
  }
  else
  {
    std::cerr << "ERROR: Illegal surface mesher specified"
                 "in chiSurfaceMesherCreate" << std::endl;
    exit(EXIT_FAILURE);
  }

  cur_hndlr.surface_mesher = new_mesher;

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "chiSurfaceMesherCreate: Surface remesher created."
    << std::endl;

  return 0;
}