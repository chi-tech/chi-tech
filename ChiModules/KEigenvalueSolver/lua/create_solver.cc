#include <ChiLua/chi_lua.h>

#include "../k_eigenvalue_solver.h"

#include <ChiPhysics//chi_physics.h>

extern ChiPhysics& chi_physics_handler;

#include "chi_log.h"

extern ChiLog& chi_log;

using namespace LinearBoltzmann;

//###################################################################
/**Creates a k-eigenvalue solver.*/
int chiKEigenvalueLBSCreateSolver(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  chi_log.Log(LOG_ALLVERBOSE_1) << "Creating K-Eigenvalue solver.";

  std::string solver_name = "KEigenvalueSolver";
  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto solver = new KEigenvalue::Solver(solver_name);

  chi_physics_handler.solver_stack.push_back(solver);

  lua_pushinteger(L,
        static_cast<lua_Integer>(chi_physics_handler.solver_stack.size() - 1));
  return 1;
}




