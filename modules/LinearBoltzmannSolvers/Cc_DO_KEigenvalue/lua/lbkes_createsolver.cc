#include "ChiLua/chi_lua.h"

#include "../lbkes_k_eigenvalue_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

using namespace lbs;

//###################################################################
/**Create the solver.
\param SolverName string Optional. A string name to use for the solver.
                         Default: `"KEigenvalueSolver"`.
 *
\return Handle A handle to the created solver.*/
int chiLBKESCreateSolver(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  chi::log.LogAllVerbose1() << "Creating k-eigenvalue solver.";

  std::string solver_name = "KEigenvalueSolver";
  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto solver = std::make_shared<DiscOrdKEigenvalueSolver>(solver_name);

  chi::object_stack.push_back(solver);

  auto n = static_cast<lua_Integer>(chi::object_stack.size() - 1);
  lua_pushinteger(L, n);
  return 1;
}
