#include <CHI_LUA/chi_lua.h>
#include <CHI_MODULES/CHI_DIFFUSION/Solver/diffusion_solver.h>
#include <CHI_PHYSICS/chi_physics.h>

extern CHI_PHYSICS chi_physics_handler;

//#############################################################################
/** Initialize the Diffusion solver.
 *
\param SolverHandle int Handle to an existing diffusion solver.

\return Success bool Returns if initialization failed.
\ingroup LuaDiffusion
\author Jan*/
/* Initialize the Diffusion solver.
 *
\param SolverHandle int Handle to an existing diffusion solver.

\return Success bool Returns if initialization failed.
\ingroup LuaDiffusion
\author Jan*/
int chiDiffusionExecute(lua_State *L)
{
  int solver_index = lua_tonumber(L,1);
  chi_diffusion::Solver* solver;

  try{
    solver = (chi_diffusion::Solver*)chi_physics_handler.solver_stack.at(solver_index);
  }
  catch(std::out_of_range o){
    std::cerr << "ERROR: Invalid solver handle." << std::endl;
    exit(EXIT_FAILURE);
  }

  solver->ExecuteS();


  return 0;
}