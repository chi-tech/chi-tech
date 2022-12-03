//
// Created by Ragusa, Jean C on 12/2/22.
//

#include "ChiConsole/chi_console.h"
#include "ChiLua/chi_lua.h"
#include "cfem_diffusion_solver.h"

double cfem_diffusion::Solver::CallLua_iXYZFunction(
  lua_State* L,
  const std::string& lua_func_name,
  const int imat,
  const chi_mesh::Vector3& xyz)
{
  //============= Load lua function
  lua_getglobal(L, lua_func_name.c_str());

  //============= Error check lua function
  if (not lua_isfunction(L, -1))
    throw std::logic_error("CallLua_iXYZFunction attempted to access lua-function, " +
                           lua_func_name + ", but it seems the function"
                                           " could not be retrieved.");

  //============= Push arguments
  lua_pushinteger(L, imat);
  lua_pushnumber(L, xyz.x);
  lua_pushnumber(L, xyz.y);
  lua_pushnumber(L, xyz.z);

  //============= Call lua function
  //4 arguments, 1 result (double), 0=original error object
  double lua_return = 0.0;
  if (lua_pcall(L,4,1,0) == 0)
  {
    LuaCheckNumberValue("CallLua_iXYZFunction", L, -1);
    lua_return = lua_tonumber(L,-1);
  }
  else
    throw std::logic_error("CallLua_iXYZFunction attempted to call lua-function, " +
                           lua_func_name + ", but the call failed." +
                           xyz.PrintStr());

  lua_pop(L,1); //pop the double, or error code

  return lua_return;
};