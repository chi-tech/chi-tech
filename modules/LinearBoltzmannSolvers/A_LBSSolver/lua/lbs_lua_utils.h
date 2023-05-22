#ifndef CHITECH_LBS_LUA_UTILS_H
#define CHITECH_LBS_LUA_UTILS_H

#include "chi_lua.h"
#include "ChiParameters/input_parameters.h"

namespace lbs
{
class LBSSolver;
}

namespace lbs::common_lua_utils
{

//void SetBoundaryOptions(LBSSolver& lbs_solver,
//                        const chi_objects::InputParameters& params);

chi_objects::InputParameters GetSyntax_SetOptions();
chi_objects::ParameterBlock
SetOptions(const chi_objects::InputParameters& params);

int chiLBSSetOptions(lua_State* L);
int chiLBSSetPhiFromFieldFunction(lua_State* L);
void RegisterLuaEntities(lua_State* L);
} // namespace lbs::common_lua_utils

#endif // CHITECH_LBS_LUA_UTILS_H
