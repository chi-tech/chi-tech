#include "ChiLua/chi_lua.h"
#include "../lbs_linear_boltzmann_solver.h"

#include "chi_runtime.h"

#include "chi_runtime.h"
#include "chi_log.h"
;

//###################################################################
/**Creates a Neutral Particle Transport solver.

\param SolverName string String name of the solver.

\return SolverHandle int Handle to the solver created.

\code
phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)
--
-- Add Groupset construction here
--
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)
--
chiLBSInitialize(phys1)
chiLBSExecute(phys1)
--
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)
\endcode


\ingroup LuaLBS
 */
int chiLBSCreateSolver(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  chi::log.LogAllVerbose1() << "Creating Linear Boltzman solver";

  std::string solver_name = "LBSolver";
  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto new_solver = std::make_shared<lbs::SteadySolver>(solver_name);

  chi::solver_stack.push_back(new_solver);

  lua_pushinteger(L,
      static_cast<lua_Integer>(chi::solver_stack.size()-1));
  return 1;
}