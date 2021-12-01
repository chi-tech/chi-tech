#include "ChiLua/chi_lua.h"

#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"

/** \defgroup LuaUnpartitionedMesh Unpartitioned Mesh-Reader
 * \ingroup LuaMesh
 */

//###################################################################
/**Creates an empty unpartitioned mesh. An empty unpartitioned mesh
 * is meant to be manipulated with calls to chiUnpartitionedMeshUploadVertex()
 * and chiUnpartitionedMeshUploadCell(). It essentially supports building a mesh
 * manually.
 *
##_

###Example
Example usage
\code
umesh = chiCreateEmptyUnpartitionedMesh()
\endcode

\ingroup LuaUnpartitionedMesh*/
int chiCreateEmptyUnpartitionedMesh(lua_State* L)
{
  const std::string func_name = __FUNCTION__;

  auto handler = chi_mesh::GetCurrentHandler();
  handler->unpartitionedmesh_stack.push_back(
    new chi_mesh::UnpartitionedMesh());

  lua_pushnumber(L,
    static_cast<lua_Number>(handler->unpartitionedmesh_stack.size()-1));

  return 1;
}

//###################################################################
/**Destroy an unpartitioned mesh. This routine should be called for
 * memory sensitive simulations because each process will have a full
 * copy of this data.
 *
\param handle int Handle to mesh.

##_

###Example
Example usage
\code
chiDestroyUnpartitionedMesh(umesh)
\endcode

\ingroup LuaUnpartitionedMesh*/
int chiDestroyUnpartitionedMesh(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);

  const int handle = lua_tointeger(L, 1);

  auto handler = chi_mesh::GetCurrentHandler();

  chi_mesh::UnpartitionedMesh* mesh_ptr;

  try{
    mesh_ptr = handler->unpartitionedmesh_stack.at(handle);
    if (mesh_ptr == nullptr) throw std::out_of_range("");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(func_name + ": Invalid mesh-handle (" +
                           std::to_string(handle) + ").");
  }

  delete mesh_ptr;
  handler->unpartitionedmesh_stack[handle] = nullptr;

  return 0;
}

//###################################################################
/**Creates an unpartitioned mesh from VTK Unstructured mesh files.

\param file_name char Filename of the .vtu file.
\param field char Name of the cell data field from which to read
                  material and boundary identifiers (optional).

\ingroup LuaUnpartitionedMesh

##_

### Example
An example mesh creation below:
\code
chiMeshHandlerCreate()

umesh = chiUnpartitionedMeshFromVTU("ZMeshTest_0.vtu")

region1 = chiRegionCreate()
chiRegionAddEmptyBoundary(region1)

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()
\endcode


\return A handle to the newly created UnpartitionedMesh*/
int chiUnpartitionedMeshFromVTU(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name,1,num_args);

  LuaCheckNilValue(func_name,L,1);
  if (num_args >= 2) LuaCheckNilValue(func_name,L,2);

  const char* temp = lua_tostring(L,1);
  const char* field = "";
  if (num_args >= 2) field = lua_tostring(L,2);
  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);
  options.material_id_fieldname = field;
  options.boundary_id_fieldname = field;

  new_object->ReadFromVTU(options);

  auto handler = chi_mesh::GetCurrentHandler();
  handler->unpartitionedmesh_stack.push_back(new_object);

  lua_pushnumber(L,
    static_cast<lua_Number>(handler->unpartitionedmesh_stack.size()-1));

  return 1;
}

//###################################################################
/**Creates an unpartitioned mesh from starccm+ exported
Ensight Gold mesh files.

\param file_name char Filename of the .case file.
\param scale float Scale to apply to the mesh

\ingroup LuaUnpartitionedMesh

##_

### Example
An example mesh creation below:
\code
chiMeshHandlerCreate()

chiUnpartitionedMeshFromEnsightGold("ChiResources/TestObjects/Sphere.case")

region1 = chiRegionCreate()
chiRegionAddEmptyBoundary(region1)

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()
\endcode

\return A handle to the newly created UnpartitionedMesh*/
int chiUnpartitionedMeshFromEnsightGold(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args <1)
    LuaPostArgAmountError(func_name,1,num_args);

  LuaCheckNilValue(func_name,L,1);
  if (num_args >= 2) LuaCheckNilValue(func_name,L,2);

  const char* temp = lua_tostring(L,1);
  double scale = 1.0;
  if (num_args >= 2) scale = lua_tonumber(L,2);
  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);
  options.scale = scale;

  new_object->ReadFromEnsightGold(options);

  auto handler = chi_mesh::GetCurrentHandler();
  handler->unpartitionedmesh_stack.push_back(new_object);

  lua_pushnumber(L,
    static_cast<lua_Number>(handler->unpartitionedmesh_stack.size()-1));

  return 1;
}

//###################################################################
/**Creates an unpartitioned mesh from a wavefront .obj file.

\param file_name char Filename of the .case file.

\ingroup LuaUnpartitionedMesh

##_

### Example
An example mesh creation below:
\code
chiMeshHandlerCreate()

chiUnpartitionedMeshFromWavefrontOBJ("ChiResources/TestObjects/TriangleMesh2x2.obj")

region1 = chiRegionCreate()
chiRegionAddEmptyBoundary(region1)

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()
\endcode

\return A handle to the newly created UnpartitionedMesh*/
int chiUnpartitionedMeshFromWavefrontOBJ(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args <1)
    LuaPostArgAmountError(func_name,1,num_args);

  LuaCheckNilValue(func_name,L,1);

  const char* temp = lua_tostring(L,1);

  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);

  new_object->ReadFromWavefrontOBJ(options);

  auto handler = chi_mesh::GetCurrentHandler();
  handler->unpartitionedmesh_stack.push_back(new_object);

  lua_pushnumber(L,
    static_cast<lua_Number>(handler->unpartitionedmesh_stack.size()-1));

  return 1;
}

//###################################################################
/**Creates an unpartitioned mesh from a .msh file.

\param file_name char Filename of the .msh file.

\ingroup LuaUnpartitionedMesh

##_

### Example
An example mesh creation below:
\code
chiMeshHandlerCreate()

chiUnpartitionedMeshFromMshFormat("File.msh")

region1 = chiRegionCreate()
chiRegionAddEmptyBoundary(region1)

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()
\endcode


\return A handle to the newly created UnpartitionedMesh*/
int chiUnpartitionedMeshFromMshFormat(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args <1)
    LuaPostArgAmountError(func_name,1,num_args);

  LuaCheckNilValue(func_name,L,1);

  const char* temp = lua_tostring(L,1);

  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);

  new_object->ReadFromMsh(options);

  auto handler = chi_mesh::GetCurrentHandler();
  handler->unpartitionedmesh_stack.push_back(new_object);

  lua_pushnumber(L,
    static_cast<lua_Number>(handler->unpartitionedmesh_stack.size()-1));

  return 1;
}
