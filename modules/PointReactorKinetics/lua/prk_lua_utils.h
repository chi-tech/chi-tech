#ifndef CHITECH_POINT_REACTOR_KINETICS_LUA_UTILS_H
#define CHITECH_POINT_REACTOR_KINETICS_LUA_UTILS_H

#include "chi_lua.h"

namespace prk::lua_utils
{
int chiPRKGetParam(lua_State* L);
int chiPRKSetParam(lua_State* L);

chi_objects::InputParameters SetParamSyntax();
chi_objects::ParameterBlock
SetParam(const chi_objects::InputParameters& params);

chi_objects::InputParameters GetParamSyntax();
chi_objects::ParameterBlock
GetParam(const chi_objects::InputParameters& params);
} // namespace prk::lua_utils

#endif