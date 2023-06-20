#ifndef CHITECH_LUA_MESHMODIFIERS_H
#define CHITECH_LUA_MESHMODIFIERS_H

#include "ChiParameters/input_parameters.h"

namespace chi_mesh::lua_utils
{
chi_objects::InputParameters MeshModifiersApply_Syntax();
chi_objects::ParameterBlock
MeshModifiersApply(const chi_objects::InputParameters& params);
} // namespace chi_mesh::lua_utils

#endif // CHITECH_LUA_MESHMODIFIERS_H
