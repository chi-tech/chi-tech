#ifndef CHITECH_UNPARTITION_MESH_LUA_UTILS_H
#define CHITECH_UNPARTITION_MESH_LUA_UTILS_H

#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

#include "chi_lua.h"

namespace chi_mesh
{
namespace unpartition_mesh_lua_utils
{
  chi_mesh::UnpartitionedMesh&
    GetUnpartitionedMeshByHandle(int handle,
                                 const std::string& calling_function_name);
}//namespace unpartition_mesh_lua_utils
}//namespace chi_mesh

//create.cc
int chiCreateEmptyUnpartitionedMesh(lua_State* L);
int chiDestroyUnpartitionedMesh(lua_State* L);

int chiUnpartitionedMeshFromVTU(lua_State* L);
int chiUnpartitionedMeshFromEnsightGold(lua_State* L);
int chiUnpartitionedMeshFromWavefrontOBJ(lua_State* L);
int chiUnpartitionedMeshFromMshFormat(lua_State* L);

//basic_operations.cc
int chiUnpartitionedMeshUploadVertex(lua_State* L);
int chiUnpartitionedMeshUploadCell(lua_State* L);
int chiUnpartitionedMeshFinalizeEmpty(lua_State* L);

#endif //CHITECH_UNPARTITION_MESH_LUA_UTILS_H
