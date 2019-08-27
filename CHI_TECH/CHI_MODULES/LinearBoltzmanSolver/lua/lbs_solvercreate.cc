#include "../../../CHI_LUA/chi_lua.h"
#include "CHI_MODULES/LinearBoltzmanSolver/lbs_linear_boltzman_solver.h"
#include "../../../CHI_PHYSICS/chi_physics.h"
#include <chi_log.h>

extern CHI_PHYSICS chi_physics_handler;
extern CHI_LOG chi_log;

/** \defgroup LuaNPT Neutral Particle Transport
 * \ingroup LuaModules*/

//###################################################################
/**Creates a Neutral Particle Transport solver.

\return SolverHandle int Handle to the solver created.

\code

\endcode


\ingroup LuaNPT
 */
int chiLBSransportCreateSolver(lua_State *L)
{
  chi_log.Log(LOG_ALLVERBOSE_1)
  << "Creating NPTransport solver";
  LinearBoltzmanSolver* new_solver = new LinearBoltzmanSolver;

  chi_physics_handler.solver_stack.push_back(new_solver);

  lua_pushnumber(L,chi_physics_handler.solver_stack.size()-1);
  return 1;
}