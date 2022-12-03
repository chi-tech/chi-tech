#include "ds_lua_utils.h"

#include "chi_runtime.h"

cfem_diffusion::Solver& cfem_diffusion::cfem_diffusion_lua_utils::
  GetSolverByHandle(int handle, const std::string& calling_function_name)
{
  std::shared_ptr<cfem_diffusion::Solver> ds_solver;
  try{

    ds_solver = std::dynamic_pointer_cast<cfem_diffusion::Solver>(
      chi::solver_stack.at(handle));

    if (not ds_solver)
      throw std::logic_error(calling_function_name +
      ": Invalid solver at given handle (" +
      std::to_string(handle) + "). "
      "The solver is not of type cfem_diffusion::Solver.");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(calling_function_name + ": Invalid solver-handle (" +
                           std::to_string(handle) + ").");
  }

  return *ds_solver;
}

#define LUA_FMACRO1(x) lua_register(L, #x, x)
#define LUA_CMACRO1(x,y) \
        lua_pushnumber(L, y); \
        lua_setglobal(L, #x)

void cfem_diffusion::cfem_diffusion_lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiCFEMDiffusionSolverCreate);
  LUA_FMACRO1(chiCFEMDiffusionSetBCProperty());

  LUA_CMACRO1(MAX_ITERATIONS, 1);
  LUA_CMACRO1(TOLERANCE     , 2);
}