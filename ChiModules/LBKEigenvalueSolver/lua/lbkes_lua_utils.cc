#include "lbkes_lua_utils.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

lbs::KEigenvalueSolver* lbs::k_eigenvalue_lua_utils::
  GetSolverByHandle(int handle, const std::string& calling_function_name)
{
  lbs::KEigenvalueSolver* lbkes_solver;
  try{

    lbkes_solver = dynamic_cast<lbs::KEigenvalueSolver*>(
      chi_physics_handler.solver_stack.at(handle));

    if (not lbkes_solver)
      throw std::logic_error(calling_function_name +
      ": Invalid solver at given handle (" +
      std::to_string(handle) + "). "
      "The solver is not of type LinearBoltzmann::KEigenvalueSolver.");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(calling_function_name + ": Invalid solver-handle (" +
                           std::to_string(handle) + ").");
  }

  return lbkes_solver;
}

#define LUA_FMACRO1(x) lua_register(L, #x, x)
#define LUA_CMACRO1(x,y) \
        lua_pushnumber(L, y); \
        lua_setglobal(L, #x)

void lbs::k_eigenvalue_lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiLBKESCreateSolver);
  LUA_FMACRO1(chiLBKESInitialize);
  LUA_FMACRO1(chiLBKESExecute);
  LUA_FMACRO1(chiLBKESSetProperty);

  LUA_CMACRO1(MAX_ITERATIONS, 1);
  LUA_CMACRO1(TOLERANCE     , 2);
}