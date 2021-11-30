#ifndef CHITECH_UNPARTITION_MESH_LUA_UTILS_H
#define CHITECH_UNPARTITION_MESH_LUA_UTILS_H

#include "ChiMesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

namespace chi_mesh
{
namespace unpartition_mesh_lua_utils
{
  chi_mesh::UnpartitionedMesh&
    GetUnpartitionedMeshByHandle(int handle,
                                 const std::string& calling_function_name);
}//namespace unpartition_mesh_lua_utils
}//namespace chi_mesh

#endif //CHITECH_UNPARTITION_MESH_LUA_UTILS_H
