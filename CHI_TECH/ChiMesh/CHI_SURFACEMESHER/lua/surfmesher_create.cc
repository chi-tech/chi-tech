#include "../../../CHI_LUA/chi_lua.h"
#include <iostream>
#include "../Predefined/surfmesher_predefined.h"
#include "../Delaunay/delaunay_mesher.h"
#include "../Triangle/triangle_mesher.h"

#include "../../CHI_MESHHANDLER/chi_meshhandler.h"

/** \defgroup LuaSurfaceMesher Surface Re-meshers
 * \ingroup LuaMesh
 *
 * chi_mesh::SurfaceMesherDelaunay
*/

#include <chi_log.h>

extern ChiLog chi_log;

//#############################################################################
/** Creates a new surface mesher remeshing.
 *
\param Type int Surface Remesher type.

Remesher types:\n
 SURFACEMESHER_PREDEFINED = No remeshing is performed.\n
 SURFACEMESHER_DELAUNAY   = Delaunay surface remesher.\n
 SURFACEMESHER_TRIANGLE   = Triangle surface remesher.

\ingroup LuaSurfaceMesher
\author Jan*/
int chiSurfaceMesherCreate(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int type = lua_tonumber(L,1);

  chi_mesh::SurfaceMesher* new_mesher;
  if (type==1)
  {
    new_mesher = new chi_mesh::SurfaceMesherPredefined;
  }
  else if (type==2)
  {
    new_mesher = new chi_mesh::SurfaceMesherDelaunay;
  }
  else if (type==3)
  {
    new_mesher = new chi_mesh::SurfaceMesherTriangle;
  } else
  {
    std::cerr << "ERROR: Illegal surface mesher specified"
                 "in chiSurfaceMesherCreate" << std::endl;
    exit(EXIT_FAILURE);
  }

  cur_hndlr->surface_mesher = new_mesher;

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "chiSurfaceMesherCreate: Surface remesher created."
    << std::endl;

  return 0;
}