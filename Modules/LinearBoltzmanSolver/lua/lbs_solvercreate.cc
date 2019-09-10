#include "ChiLua/chi_lua.h"
#include "../lbs_linear_boltzman_solver.h"
#include "ChiPhysics/chi_physics.h"
#include <chi_log.h>

extern ChiPhysics chi_physics_handler;
extern ChiLog chi_log;

/** \defgroup LuaNPT Linear Boltzman Solver
 * \ingroup LuaModules*/

//###################################################################
/**Creates a Neutral Particle Transport solver.

\return SolverHandle int Handle to the solver created.

\code
phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)
--
-- Add Groupset construction here
--
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)
--
chiLBSInitialize(phys1)
chiLBSExecute(phys1)
--
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)
\endcode


\ingroup LuaNPT
 */
int chiLBSCreateSolver(lua_State *L)
{
  chi_log.Log(LOG_ALLVERBOSE_1)
  << "Creating Linear Boltzman solver";
  LinearBoltzmanSolver* new_solver = new LinearBoltzmanSolver;

  chi_physics_handler.solver_stack.push_back(new_solver);

  lua_pushnumber(L,chi_physics_handler.solver_stack.size()-1);
  return 1;
}