#include "chi_lua.h"

#include "ChiConsole/chi_console.h"

#include "../functionxyz_dimA_to_dimB.h"

#include "chi_log_exceptions.h"

namespace chi_math::functions::lua_utils
{

int chiFunctionXYZDimAToDimBEvaluate(lua_State* L);

RegisterLuaFunctionAsIs(chiFunctionXYZDimAToDimBEvaluate);

/**Evaluates a function of base type `FunctionXYZDimAToDimB`.
\param handle int. Handle to the function to evaluate.
\param x double. The x position.
\param y double. The y position.
\param z double. The z position.
\param phi_vals table. This array of doubles are the input values. They could
                       coupled variables.

\return Table A table of output values.*/
int chiFunctionXYZDimAToDimBEvaluate(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 5) LuaPostArgAmountError(fname, 5, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNumberValue(fname, L, 2);
  LuaCheckNumberValue(fname, L, 3);
  LuaCheckNumberValue(fname, L, 4);
  LuaCheckTableValue(fname, L, 5);

  const size_t handle = lua_tointeger(L, 1);
  const double x = lua_tonumber(L, 2);
  const double y = lua_tonumber(L, 3);
  const double z = lua_tonumber(L, 4);

  chi_objects::ParameterBlock table =
    chi_lua::TableParserAsParameterBlock::ParseTable(L, 5);

  const auto array_vals = table.GetVectorValue<double>();

  if (table.Type() != chi_objects::ParameterBlockType::ARRAY)
    ChiInvalidArgument("Parameter 5 is required to be an array. This means"
                       "the keys cannot have names.");

  const auto& function =
    chi::GetStackItem<chi_math::functions::FunctionXYZDimAToDimB>(
      chi::object_stack, handle, fname);

  const std::vector<double> values = function.Evaluate(x, y, z, array_vals);

  if (values.size() == 1)
  {
    lua_pushnumber(L, values.front());
    return 1;
  }
  // else

  lua_newtable(L);
  for (size_t k = 0; k < values.size(); ++k)
  {
    lua_pushinteger(L, static_cast<lua_Integer>(k)+1);
    lua_pushnumber(L, values[k]);
    lua_settable(L, -3);
  }
  return 1;
}

} // namespace chi_math::functions::lua_utils