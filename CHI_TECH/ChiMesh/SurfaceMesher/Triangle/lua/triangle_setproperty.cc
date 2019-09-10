#include "../../../../ChiLua/chi_lua.h"
#include <iostream>
#include "../triangle_mesher.h"

#include "../../../MeshHandler/chi_meshhandler.h"
#include <chi_log.h>

extern ChiLog chi_log;

//#############################################################################
/** Sets a property of a Triangle surface mesher.

\param PropertyNumber int Handle of the property to be set.
\param PropertyValue varying Value of the property.

Properties:\n
 MAX_AREA = Area constraint.\n
 PARTITION_X   = Number of partitions in X.\n
 PARTITION_Y   = Number of partitions in Y.\n
 CUT_X = Adds a cut at the given x-value.\n
 CUT_Y = Adds a cut at the given y-value.

\ingroup LuaSurfaceMesher
\author Jan*/
int chiSurfaceMesherTriangleSetProperty(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  if (dynamic_cast<chi_mesh::SurfaceMesherTriangle*>
                  (cur_hndlr->surface_mesher) == NULL )
  {
    std::cerr << "ERROR! Current surface mesher is not"
                 "of type SURFACEMESHER_TRIANGLE.\n";
    exit(EXIT_FAILURE);
  }

  chi_mesh::SurfaceMesherTriangle* tri_mesher =
    (chi_mesh::SurfaceMesherTriangle*)(cur_hndlr->surface_mesher);

  //================================================== Get property number
  int property_num = lua_tonumber(L,1);

  //================================================== Area constraint
  if (property_num == 1)   //MAX_AREA
  {
    double constraint = lua_tonumber(L,2);
    tri_mesher->area_constraint = constraint;
    tri_mesher->get_auto_min = false;
    chi_log.Log(LOG_0VERBOSE_2) << "SurfaceMesherTriangle: "
                                 "Area constraint set to "
                              << constraint;
  }
  //================================================== Partitioning
  if (property_num == 2)   //PARTITION_X
  {
    int num = lua_tonumber(L,2);
    tri_mesher->partitioning_x = num;
  }
  if (property_num == 3)   //PARTITION_Y
  {
    int num = lua_tonumber(L,2);
    tri_mesher->partitioning_y = num;
  }
  if (property_num == 4)   //CUT_X
  {
    double cut = lua_tonumber(L,2);
    tri_mesher->xcuts.push_back(cut);
  }
  if (property_num == 5)   //CUT_Y
  {
    double cut = lua_tonumber(L,2);
    tri_mesher->ycuts.push_back(cut);
  }



  return 0;
}