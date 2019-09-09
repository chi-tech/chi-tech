#include"ChiLua/chi_lua.h"

#include"../Solver/diffusion_solver.h"

#include"ChiPhysics/chi_physics.h"
#include <chi_log.h>

extern ChiPhysics chi_physics_handler;
extern ChiLog chi_log;


//#############################################################################
/** Creates a Diffusion solver.

\return Handle int Handle to the created solver.
\ingroup LuaDiffusion
\author Jan*/
int chiDiffusionCreateSolver(lua_State *L)
{
  chi_diffusion::Solver* new_solver = new chi_diffusion::Solver;

  chi_physics_handler.solver_stack.push_back(new_solver);

  lua_pushnumber(L,chi_physics_handler.solver_stack.size()-1);

  chi_log.Log(LOG_ALLVERBOSE_1)
    << "chiDiffusionCreateSolver: Diffusion solver created"
    << std::endl;
  return 1;
}