#include "unpartition_mesh_lua_utils.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

namespace chi_mesh
{
namespace unpartition_mesh_lua_utils
{
chi_mesh::UnpartitionedMesh&
  GetUnpartitionedMeshByHandle(int handle,
                               const std::string& calling_function_name)
{
  const std::string func_name = __FUNCTION__;

  auto handler = chi_mesh::GetCurrentHandler();

  chi_mesh::UnpartitionedMesh* mesh_ptr;

  try{
    mesh_ptr = handler->unpartitionedmesh_stack.at(handle);
    if (mesh_ptr == nullptr) throw std::out_of_range("");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(calling_function_name + ": Invalid mesh-handle (" +
                           std::to_string(handle) + ").");
  }

  return *mesh_ptr;
}
}//namespace unpartition_mesh_lua_utils
}//namespace chi_mesh
