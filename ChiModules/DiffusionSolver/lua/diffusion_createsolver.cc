#include"ChiLua/chi_lua.h"

#include"../Solver/diffusion_solver.h"

#include"ChiPhysics/chi_physics.h"
#include <chi_log.h>

extern ChiPhysics&  chi_physics_handler;
extern ChiLog& chi_log;


//#############################################################################
/** Creates a Diffusion solver.

\return Handle int Handle to the created solver.
\ingroup LuaDiffusion
\author Jan*/
int chiDiffusionCreateSolver(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  std::string solver_name = "DiffusionSolver";

  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto new_solver = new chi_diffusion::Solver(solver_name);

  chi_physics_handler.solver_stack.push_back(new_solver);

  lua_pushinteger(L,
      static_cast<lua_Integer>(chi_physics_handler.solver_stack.size()-1));

  chi_log.Log(LOG_ALLVERBOSE_1)
    << "chiDiffusionCreateSolver: Diffusion solver created"
    << std::endl;
  return 1;
}