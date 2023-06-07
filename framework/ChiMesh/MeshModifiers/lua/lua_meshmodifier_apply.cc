#include "lua_meshmodifiers.h"

#include "chi_runtime.h"

#include "ChiMesh/MeshModifiers/MeshModifier.h"

#include "ChiConsole/chi_console.h"

namespace chi_mesh::lua_utils
{

RegisterWrapperFunction(/*namespace_in_lua=*/chi_mesh,
                        /*name_in_lua=*/MeshModifiersApply,
                        /*syntax_function=*/MeshModifiersApply_Syntax,
                        /*actual_function=*/MeshModifiersApply);

chi_objects::InputParameters MeshModifiersApply_Syntax()
{
  chi_objects::InputParameters params;

  params.SetGeneralDescription(
    "\\defgroup chi_mesh__MeshModifiersApply chi_mesh.MeshModifiersApply \n"
    "\\ingroup DocMeshModifiers\n"
    "Lua wrapper function for applying mesh modifiers");

  params.AddRequiredParameterArray(
    "arg0", "A list of handles to the modifiers to apply.");

  return params;
}

chi_objects::ParameterBlock
MeshModifiersApply(const chi_objects::InputParameters& params)
{
  const std::string fname = __FUNCTION__;
  const std::vector<size_t> handles =
    params.GetParamVectorValue<size_t>("arg0");

  for (const size_t handle : handles)
  {
    auto& modifier =
      chi::GetStackItem<MeshModifier>(chi::object_stack, handle, fname);

    modifier.Apply();
  }

  return chi_objects::ParameterBlock(); // Return empty param block
}

} // namespace chi_mesh::lua_utils