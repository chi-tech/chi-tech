#include "ChiLua/chi_lua.h"

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "chi_log.h"

extern ChiLog&  chi_log;

/** \defgroup LuaMeshMacros Mesh Macros
 * These functions are considered "macros" because they encapsulate
 * functionality available from lower level function calls.
 *
 * \ingroup LuaMesh
*/

//###################################################################
/** Creates a 1D Mesh from an array of 1D vertices.

\param x_nodes array_float An Array of floating point numbers denoting
                           1D nodes along x-axis.

\return Two handles: line-mesh, region

\ingroup LuaMeshMacros

##_

### Example
An example 1D mesh creation below:
\code
chiMeshHandlerCreate()
nodes={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
chiMeshCreate1DSlabMesh(nodes)
chiVolumeMesherExecute();
\endcode

 \author Nak*/
int chiMeshCreate1DSlabMesh(lua_State* L)
{
  //=================================== Check argc
  const char func_name[] = "chiMeshCreate1DSlabMesh";
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(func_name,1,num_args);

  //=================================== Check args table
  if (not lua_istable(L,1))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": First argument found to not be an array.";
    exit(EXIT_FAILURE);
  }

  //=================================== Decl vars
  int table_index=0;
  int N=0;
  std::vector<std::vector<double>> array(3);

  //=================================== Get first array
  table_index = 1;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  //=================================== Create mesh
  chi_mesh::Create1DSlabMesh(array[0]);

  //=================================== Push handles
  auto handler = chi_mesh::GetCurrentHandler();
  lua_pushnumber(L,handler->linemesh_stack.size()-1);
  lua_pushnumber(L,handler->region_stack.size()-1);

  return 2;
}

//###################################################################
/** Creates a 1D Mesh from an array of 1D vertices.

\param x_nodes array_float An Array of floating point numbers denoting
                           1D nodes along x-axis.

\return Two handles: line-mesh, region

\ingroup LuaMeshMacros

##_

### Example
An example 1D mesh creation below:
\code
chiMeshHandlerCreate()
nodes={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
chiMeshCreateUnpartitioned1DOrthoMesh(nodes)
chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
chiVolumeMesherExecute();
\endcode

 \author Nak*/
int chiMeshCreateUnpartitioned1DOrthoMesh(lua_State* L)
{
  //=================================== Check argc
  const char func_name[] = "chiMeshCreateUnpartitioned1DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(func_name,1,num_args);

  //=================================== Check args table
  if (not lua_istable(L,1))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": First argument found to not be an array.";
    exit(EXIT_FAILURE);
  }

  //=================================== Decl vars
  int table_index=0;
  int N=0;
  std::vector<std::vector<double>> array(3);

  //=================================== Get first array
  table_index = 1;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  //=================================== Create mesh
  chi_mesh::CreateUnpartitioned1DOrthoMesh(array[0]);

  //=================================== Push handles
  auto handler = chi_mesh::GetCurrentHandler();
  lua_pushnumber(L,handler->unpartitionedmesh_stack.size()-1);
  lua_pushnumber(L,handler->region_stack.size()-1);

  return 2;
}

//###################################################################
/** Creates a 2D Orthogonal Mesh from arrays of 1D vertices.

\param x_nodes array_float An Array of floating point numbers denoting
                           1D nodes along x-axis.
\param y_nodes array_float An Array of floating point numbers denoting
                           1D nodes along y-axis.

\return Two handles: surface-mesh, region

\ingroup LuaMeshMacros

##_

### Example
An example 2D mesh creation below:
\code
chiMeshHandlerCreate()
nodesx={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesy={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
chiMeshCreate2DOrthoMesh(nodesx,nodesy)
chiVolumeMesherExecute();
\endcode

 \author Nak*/
int chiMeshCreate2DOrthoMesh(lua_State* L)
{
  //=================================== Check argc
  const char func_name[] = "chiMeshCreate2DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(func_name,2,num_args);

  //=================================== Check args table
  if (not lua_istable(L,1))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": First argument found to not be an array.";
    exit(EXIT_FAILURE);
  }
  if (not lua_istable(L,2))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": Second argument found to not be an array.";
    exit(EXIT_FAILURE);
  }

  //=================================== Decl vars
  int table_index=0;
  int N=0;
  std::vector<std::vector<double>> array(3);

  //=================================== Get first array
  table_index = 1;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  //=================================== Get second array
  table_index = 2;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  //=================================== Create mesh
  chi_mesh::Create2DOrthoMesh(array[0],array[1]);

  //=================================== Push handles
  auto handler = chi_mesh::GetCurrentHandler();
  lua_pushnumber(L,handler->surface_mesh_stack.size()-1);
  lua_pushnumber(L,handler->region_stack.size()-1);

  return 2;
}

//###################################################################
/** Creates a 2D Orthogonal Mesh from arrays of 1D vertices.

\param x_nodes array_float An Array of floating point numbers denoting
                           1D nodes along x-axis.
\param y_nodes array_float An Array of floating point numbers denoting
                           1D nodes along y-axis.

\return Two handles: surface-mesh, region

\ingroup LuaMeshMacros

##_

### Example
An example 2D mesh creation below:
\code
chiMeshHandlerCreate()
nodesx={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesy={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
chiMeshCreateUnpartitioned2DOrthoMesh(nodesx,nodesy)
chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
chiVolumeMesherExecute();
\endcode

 \author Nak*/
int chiMeshCreateUnpartitioned2DOrthoMesh(lua_State* L)
{
  //=================================== Check argc
  const char func_name[] = "chiMeshCreateUnpartitioned2DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(func_name,2,num_args);

  //=================================== Check args table
  if (not lua_istable(L,1))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": First argument found to not be an array.";
    exit(EXIT_FAILURE);
  }
  if (not lua_istable(L,2))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": Second argument found to not be an array.";
    exit(EXIT_FAILURE);
  }

  //=================================== Decl vars
  int table_index=0;
  int N=0;
  std::vector<std::vector<double>> array(3);

  //=================================== Get first array
  table_index = 1;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  //=================================== Get second array
  table_index = 2;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  //=================================== Create mesh
  chi_mesh::CreateUnpartitioned2DOrthoMesh(array[0],array[1]);

  //=================================== Push handles
  auto handler = chi_mesh::GetCurrentHandler();
  lua_pushnumber(L,handler->unpartitionedmesh_stack.size()-1);
  lua_pushnumber(L,handler->region_stack.size()-1);

  return 2;
}

//###################################################################
/** Creates a 3D Orthogonal Mesh from arrays of 1D vertices. The
 * underlying mesher is an extruder.

\param x_nodes array_float An Array of floating point numbers denoting
                           1D nodes along x-axis.
\param y_nodes array_float An Array of floating point numbers denoting
                           1D nodes along y-axis.
\param z_nodes array_float An Array of floating point numbers denoting
                           1D nodes along z-axis.

\return Two handles: surface-mesh, region

\ingroup LuaMeshMacros

##_

### Example
An example 3D mesh creation below:
\code
chiMeshHandlerCreate()
nodesx={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesy={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesz={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
chiMeshCreate3DOrthoMesh(nodesx,nodesy,nodesz)
chiVolumeMesherExecute();
\endcode

 \author Nak*/
int chiMeshCreate3DOrthoMesh(lua_State* L)
{
  //=================================== Check argc
  const char func_name[] = "chiMeshCreate3DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError(func_name,3,num_args);

  //=================================== Check args table
  if (not lua_istable(L,1))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": First argument found to not be an array.";
    exit(EXIT_FAILURE);
  }
  if (not lua_istable(L,2))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": Second argument found to not be an array.";
    exit(EXIT_FAILURE);
  }
  if (not lua_istable(L,3))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": Third argument found to not be an array.";
    exit(EXIT_FAILURE);
  }

  //=================================== Decl vars
  int table_index=0;
  int N=0;
  std::vector<std::vector<double>> array(3);

  //=================================== Get first array
  table_index = 1;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  //=================================== Get second array
  table_index = 2;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  //=================================== Get second array
  table_index = 3;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  //=================================== Create mesh
  chi_mesh::Create3DOrthoMesh(array[0],array[1],array[2]);

  //=================================== Push handles
  auto handler = chi_mesh::GetCurrentHandler();
  lua_pushnumber(L,handler->surface_mesh_stack.size()-1);
  lua_pushnumber(L,handler->region_stack.size()-1);

  return 2;
}

//###################################################################
/** Creates a 3D Orthogonal Mesh from arrays of 1D vertices. The
 * underlying mesher is an extruder.

\param x_nodes array_float An Array of floating point numbers denoting
                           1D nodes along x-axis.
\param y_nodes array_float An Array of floating point numbers denoting
                           1D nodes along y-axis.
\param z_nodes array_float An Array of floating point numbers denoting
                           1D nodes along z-axis.

\return Two handles: surface-mesh, region

\ingroup LuaMeshMacros

##_

### Example
An example 3D mesh creation below:
\code
chiMeshHandlerCreate()
nodesx={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesy={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesz={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
chiMeshCreateUnpartitioned3DOrthoMesh(nodesx,nodesy,nodesz)
chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
chiVolumeMesherExecute();
\endcode

 \author Nak*/
int chiMeshCreateUnpartitioned3DOrthoMesh(lua_State* L)
{
  //=================================== Check argc
  const char func_name[] = "chiMeshCreateUnpartitioned3DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError(func_name,3,num_args);

  //=================================== Check args table
  if (not lua_istable(L,1))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": First argument found to not be an array.";
    exit(EXIT_FAILURE);
  }
  if (not lua_istable(L,2))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": Second argument found to not be an array.";
    exit(EXIT_FAILURE);
  }
  if (not lua_istable(L,3))
  {
    chi_log.Log(LOG_ALLERROR)
      << func_name << ": Third argument found to not be an array.";
    exit(EXIT_FAILURE);
  }

  //=================================== Decl vars
  int table_index=0;
  int N=0;
  std::vector<std::vector<double>> array(3);

  //=================================== Get first array
  table_index = 1;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  //=================================== Get second array
  table_index = 2;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  //=================================== Get second array
  table_index = 3;
  N = lua_rawlen(L,table_index);
  array[table_index-1].resize(N);
  for (int k=0; k<N; k++)
  {
    lua_pushnumber(L,k+1);
    lua_gettable(L,table_index);

    array[table_index-1][k] = lua_tonumber(L,-1);
    lua_pop(L,1);
  }

  //=================================== Create mesh
  chi_mesh::CreateUnpartitioned3DOrthoMesh(array[0],array[1],array[2]);

  //=================================== Push handles
  auto handler = chi_mesh::GetCurrentHandler();
  lua_pushnumber(L,handler->unpartitionedmesh_stack.size()-1);
  lua_pushnumber(L,handler->region_stack.size()-1);

  return 2;
}