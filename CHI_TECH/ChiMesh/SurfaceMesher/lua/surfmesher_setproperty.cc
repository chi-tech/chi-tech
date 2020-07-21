#include "../../../ChiLua/chi_lua.h"
#include <iostream>
#include "../surfacemesher.h"

#include "../../MeshHandler/chi_meshhandler.h"
#include <chi_log.h>

extern ChiLog& chi_log;

//#############################################################################
/** Sets a property of a surface mesher.

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
int chiSurfaceMesherSetProperty(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  chi_mesh::SurfaceMesher* surf_mesher = cur_hndlr->surface_mesher;

  //================================================== Get property number
  int property_num = lua_tonumber(L,1);

  //================================================== Area constraint
  if (property_num == 1)   //MAX_AREA
  {
    chi_log.Log(LOG_0WARNING) << "Deprecated and removed feature"
                                 "property MAX_AREA in call"
                                 " to chiSurfaceMesherSetProperty";
  }
  //================================================== Partitioning
  if (property_num == 2)   //PARTITION_X
  {
    int num = lua_tonumber(L,2);
    surf_mesher->partitioning_x = num;

    chi_log.Log(LOG_0) << "Surface mesher partitioning x set to: " << num;
  }
  if (property_num == 3)   //PARTITION_Y
  {
    int num = lua_tonumber(L,2);
    surf_mesher->partitioning_y = num;

    chi_log.Log(LOG_0) << "Surface mesher partitioning y set to: " << num;
  }
  if (property_num == 4)   //CUT_X
  {
    double cut = lua_tonumber(L,2);
    surf_mesher->xcuts.push_back(cut);
  }
  if (property_num == 5)   //CUT_Y
  {
    double cut = lua_tonumber(L,2);
    surf_mesher->ycuts.push_back(cut);
  }



  return 0;
}