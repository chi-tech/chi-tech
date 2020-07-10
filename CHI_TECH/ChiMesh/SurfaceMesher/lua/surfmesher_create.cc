#include "../../../ChiLua/chi_lua.h"
#include <iostream>
#include "../Predefined/surfmesher_predefined.h"
#include "../PassThrough/surfmesher_passthrough.h"
#include "../Delaunay/delaunay_mesher.h"
#include "../Triangle/triangle_mesher.h"

#include "../../MeshHandler/chi_meshhandler.h"

/** \defgroup LuaSurfaceMesher Surface Re-meshers
 * \ingroup LuaMesh
 *
 * chi_mesh::SurfaceMesherDelaunay
*/

#include <chi_log.h>

extern ChiLog& chi_log;

//#############################################################################
/** Creates a surface preprocessor.
 *
\param Type int Surface Remesher type. See SurfaceMesherType.

## _

###SurfaceMesherType:\n
SurfaceMesherType.Passthrough\n
 Makes no modification to the region surfaces.\n\n

\code
chiSurfaceMesherCreate(SurfaceMesherType.Passthrough)
\endcode

SurfaceMesherType.Delaunay:\n
 Experimental. Performs a Delaunay triangulation of the region surfaces.

## _

### Legacy

 SURFACEMESHER_PREDEFINED = No remeshing is performed.\n
 SURFACEMESHER_DELAUNAY   = Delaunay surface remesher.\n
 SURFACEMESHER_TRIANGLE   = Triangle surface remesher.

\ingroup LuaSurfaceMesher
\author Jan*/
int chiSurfaceMesherCreate(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Get argument
  LuaCheckNilValue("chiSurfaceMesherCreate",L,1);
  int type = lua_tonumber(L,1);

  //============================================= Create the surface mesher
  chi_mesh::SurfaceMesher* new_mesher;
  if (type==(int)chi_mesh::SurfaceMesherType::Passthrough)
  {
    new_mesher = new chi_mesh::SurfaceMesherPredefined;
  }
  else if (type==(int)chi_mesh::SurfaceMesherType::Delaunay)
  {
    new_mesher = new chi_mesh::SurfaceMesherDelaunay;
  }
//  else if (type==3)
//  {
//    new_mesher = new chi_mesh::SurfaceMesherTriangle;
//  }
  else
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