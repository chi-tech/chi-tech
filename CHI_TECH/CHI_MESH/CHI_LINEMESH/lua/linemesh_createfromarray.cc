#include"../../../CHI_LUA/chi_lua.h"
#include<iostream>
#include "../chi_linemesh.h"
#include "../../CHI_MESHHANDLER/chi_meshhandler.h"

/** \defgroup LuaLineMesh Line Meshes
 * \ingroup LuaMesh
*/

#include <chi_log.h>

extern CHI_LOG chi_log;

//#############################################################################
/** Creates a new line mesh from a lua array.
 *
\param Table LuaTable A lua table containing all the vertices of the mesh.

## Example
\code
mesh={}
mesh[1] = {0.1,0.001}
mesh[2] = {0.2,0.001,0.003}
line_mesh = chiLineMeshCreateFromArray(mesh)
\endcode

\return Handle int Handle to the created line mesh.
\ingroup LuaLineMesh
\author Jan*/
int chiLineMeshCreateFromArray(lua_State *L)
{
  int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError("chiLineMeshCreateFromTable",1,num_args);

  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Create LineMesh
  chi_mesh::LineMesh* new_line = new chi_mesh::LineMesh;

  if (!lua_istable(L,1))
  {
    chi_log.Log(LOG_ALLERROR)
      << "In call to chiLineMeshCreateFromArray: "
      << "The arguments was detected not to be a lua table.";
    exit(EXIT_FAILURE);
  }

  //============================================= Extract values from table
  int table_len = lua_rawlen(L,1);
  std::vector<double> values_x(table_len,0.0);
  std::vector<double> values_y(table_len,0.0);
  std::vector<double> values_z(table_len,0.0);

  for (int g=0; g<table_len; g++)
  {
    lua_pushnumber(L,g+1);
    lua_gettable(L,1);

    if (lua_isnumber(L,-1))
    {
      values_x[g] = 0.0;
      values_y[g] = 0.0;
      values_z[g] = lua_tonumber(L,-1);
      lua_pop(L,1);
    }
    else if (lua_istable(L,-1))
    {
      int vec_length = lua_rawlen(L,-1);

      if (vec_length != 3)
      {
        chi_log.Log(LOG_ALLERROR)
          << "In call to chiLineMeshCreateFromArray: "
          << "Value " << g+1 << " was detected to be of length "
          << vec_length << " but a length of 3 is"
          << " required for building a vector.";
        exit(EXIT_FAILURE);
      }

      //Get x value
      lua_pushnumber(L,1);
      lua_gettable(L,2);

      values_x[g] = lua_tonumber(L,-1);
      lua_pop(L,1);

      //Get y value
      lua_pushnumber(L,2);
      lua_gettable(L,2);

      values_y[g] = lua_tonumber(L,-1);
      lua_pop(L,1);

      //Get z value
      lua_pushnumber(L,3);
      lua_gettable(L,2);

      values_z[g] = lua_tonumber(L,-1);
      lua_pop(L,1);

      lua_pop(L,1); //Pop the vector off
    }

    new_line->vertices.push_back(chi_mesh::Vector(values_x[g],
                                                  values_y[g],
                                                  values_z[g]));
  }

  //============================================= Add to handler
  cur_hndlr->linemesh_stack.push_back(new_line);

  int index = cur_hndlr->linemesh_stack.size()-1;
  lua_pushnumber(L,index);

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "chiLineMeshCreateFromTable: Created line mesh " << index << std::endl;

  return 1;
}