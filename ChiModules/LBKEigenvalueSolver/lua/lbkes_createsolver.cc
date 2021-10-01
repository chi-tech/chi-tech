#include "../lbkes_k_eigenvalue_solver.h"

#include <chi_lua.h>

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics& chi_physics_handler;

#include <chi_log.h>
extern ChiLog& chi_log;

using namespace LinearBoltzmann;

//###################################################################
/**Create the solver.*/
int chiLBKESCreateSolver(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  chi_log.Log(LOG_ALLVERBOSE_1) << "Creating k-eigenvalue solver.";

  std::string solver_name = "KEigenvalueSolver";
  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto solver = new KEigenvalueSolver(solver_name);

  chi_physics_handler.solver_stack.push_back(solver);

  auto n = static_cast<lua_Integer>(chi_physics_handler.solver_stack.size() - 1);
  lua_pushinteger(L, n);
  return 1;
}
