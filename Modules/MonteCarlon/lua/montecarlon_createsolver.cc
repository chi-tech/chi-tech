#include"ChiLua/chi_lua.h"

#include"../Solver/solver_montecarlon.h"

#include"ChiPhysics/chi_physics.h"

extern ChiPhysics chi_physics_handler;


//#############################################################################
/** Creates a MonteCarlon solver.

\return Handle int Handle to the created solver.
\ingroup LuaMonteCarlon
\author Jan*/
int chiMonteCarlonCreateSolver(lua_State *L)
{
  chi_montecarlon::Solver* new_solver = new chi_montecarlon::Solver;

  chi_physics_handler.solver_stack.push_back(new_solver);

  lua_pushnumber(L,chi_physics_handler.solver_stack.size()-1);

  return 1;
}
