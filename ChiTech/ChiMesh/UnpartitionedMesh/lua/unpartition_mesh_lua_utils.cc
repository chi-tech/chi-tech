#include "unpartition_mesh_lua_utils.h"

#include "chi_runtime.h"

namespace chi_mesh::unpartition_mesh_lua_utils
{
chi_mesh::UnpartitionedMesh&
  GetUnpartitionedMeshByHandle(int handle,
                               const std::string& calling_function_name)
{
  const std::string func_name = __FUNCTION__;

  return chi::GetStackItem<chi_mesh::UnpartitionedMesh>(
    chi::unpartitionedmesh_stack,
    handle, func_name);
}
}//namespace chi_mesh
