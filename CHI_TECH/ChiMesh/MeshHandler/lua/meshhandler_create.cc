#include"../../../ChiLua/chi_lua.h"
#include <iostream>

#include "../chi_meshhandler.h"
#include <chi_log.h>

extern ChiLog& chi_log;

extern std::vector<chi_mesh::MeshHandler*>  chi_meshhandler_stack;
extern int                                  chi_current_mesh_handler;
/** \defgroup LuaMeshHandler Mesh Handler
 * \ingroup LuaMesh
*/
//#############################################################################
/** Creates a mesh handler and sets it as "current".

\return Handle int Handle to the created mesh handler.
\ingroup LuaMeshHandler
\author Jan*/
int chiMeshHandlerCreate(lua_State *L)
{
  chi_mesh::MeshHandler* new_handler = new chi_mesh::MeshHandler;

  chi_meshhandler_stack.push_back(new_handler);

  int index = chi_meshhandler_stack.size()-1;
  lua_pushnumber(L,index);
  chi_current_mesh_handler = index;

  chi_log.Log(LOG_ALLVERBOSE_2)
  << "chiMeshHandlerCreate: Mesh Handler " << index << " created\n";

  return 1;
}


//#############################################################################
/** Sets the given mesh handler as "current".

\param HandlerHandler int Handle to the mesh handler previously created
       with a call to chiMeshHandlerCreate.

\ingroup LuaMeshHandler
\author Jan*/
int chiMeshHandlerSetCurrent(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiMeshHandlerSetCurrent",1,num_args);

  int handle = lua_tonumber(L,1);

  if ((handle < 0) or (handle >= chi_meshhandler_stack.size()))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid handle to mesh handler specified "
      << "in call to chiMeshHandlerSetCurrent";
    exit(EXIT_FAILURE);
  }

  chi_current_mesh_handler = handle;

  chi_log.Log(LOG_ALLVERBOSE_2)
    << "chiMeshHandlerSetCurrent: set to " << handle;

  return 0;
}