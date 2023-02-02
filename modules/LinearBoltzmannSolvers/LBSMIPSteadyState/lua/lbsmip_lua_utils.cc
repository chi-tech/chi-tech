#include "lbsmip_lua_utils.h"

#include "../lbsmip_steady_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace lbs::lbsmip_lua_utils
{

//###################################################################
/**Create the solver.
\param SolverName string Optional. A string name to use for the solver.
                         [Default="MIPSteadyStateSolver"].
 *
\return Handle A handle to the created solver.*/
int chiLBSMIPCreateSolver(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  chi::log.LogAllVerbose1() << "Creating MIPSteadyStateSolver solver.";

  std::string solver_name = "MIPSteadyStateSolver";
  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto solver = std::make_shared<MIPSteadyStateSolver>(solver_name);

  chi::solver_stack.push_back(solver);

  auto n = static_cast<lua_Integer>(chi::solver_stack.size() - 1);
  lua_pushinteger(L, n);
  return 1;
}

void RegisterLuaEntities(lua_State* L)
{
  lua_register(L, "chiLBSMIPCreateSolver",
               lbs::lbsmip_lua_utils::chiLBSMIPCreateSolver);
}

}//lbs::lbsmip_lua_utils