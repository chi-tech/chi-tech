#include "chi_lua.h"
#include "ChiConsole/chi_console.h"

#include "../point_reactor_kinetics.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace prk::lua_utils
{
int chiPRKGetParam(lua_State* L);
int chiPRKSetParam(lua_State* L);
ChiConsoleRegisterLuaFunction(chiPRKGetParam);
ChiConsoleRegisterLuaFunction(chiPRKSetParam);

/**Gets a parameter from the prk::TransientSolver.
*
* \param handle int Handle of the solver.
* \param param_name  string Name of the parameter to retrieve.
\return Varying*/
int chiPRKGetParam(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, num_args, 2);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);

  const int handle = lua_tointeger(L, 1);

  auto solver =
    chi::GetStackItem<TransientSolver>(chi::object_stack, handle, fname);

  const std::string param_name = lua_tostring(L, 2);

  if (param_name == "population_prev")
  {
    lua_pushnumber(L, solver.PopulationPrev());
    return 1;
  }
  else if (param_name == "population_next")
  {
    lua_pushnumber(L, solver.PopulationNext());
    return 1;
  }
  else if (param_name == "period")
  {
    lua_pushnumber(L, solver.Period());
    return 1;
  }
  else if (param_name == "time_prev")
  {
    lua_pushnumber(L, solver.TimePrev());
    return 1;
  }
  else if (param_name == "time_next")
  {
    lua_pushnumber(L, solver.TimeNext());
    return 1;
  }
  else
    throw std::invalid_argument(fname + ": Invalid parameter \"" + param_name +
                                "\".");

  return 0;
}

/**Gets a parameter from the prk::TransientSolver.
*
* \param handle int Handle of the solver.
* \param param_name  string Name of the parameter to retrieve.
* \param value Varying The value to be set to the parameter.
\return Varying*/
int chiPRKSetParam(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, num_args, 3);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  const int handle = lua_tointeger(L, 1);

  auto& solver =
    chi::GetStackItem<TransientSolver>(chi::object_stack, handle, fname);

  const std::string param_name = lua_tostring(L, 2);

  if (param_name == "rho")
  {
    LuaCheckNumberValue(
      fname + "(handle,\"rho\", : Expects a number value.", L, 3);
    const double val = lua_tonumber(L, 3);
    solver.SetRho(val);
  }
  else
    throw std::invalid_argument(fname + ": Invalid parameter \"" + param_name +
                                "\".");

  return 0;
}

} // namespace prk::lua_utils