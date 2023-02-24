#include "unpartition_mesh_lua_utils.h"

#define RegisterFunction(x) \
  lua_register(L, #x, chi_mesh::unpartition_mesh_lua_utils::x)

namespace chi_mesh::unpartition_mesh_lua_utils
{
void RegisterLuaEntities(lua_State* L)
{
  RegisterFunction(chiCreateEmptyUnpartitionedMesh);
  RegisterFunction(chiDestroyUnpartitionedMesh);

  RegisterFunction(chiUnpartitionedMeshFromVTU);
  RegisterFunction(chiUnpartitionedMeshFromPVTU);
  RegisterFunction(chiUnpartitionedMeshFromEnsightGold);
  RegisterFunction(chiUnpartitionedMeshFromWavefrontOBJ);
  RegisterFunction(chiUnpartitionedMeshFromMshFormat);
  RegisterFunction(chiUnpartitionedMeshFromExodusII);

  RegisterFunction(chiUnpartitionedMeshUploadVertex);
  RegisterFunction(chiUnpartitionedMeshUploadCell);
  RegisterFunction(chiUnpartitionedMeshFinalizeEmpty);
}
}//namespace chi_mesh
