#include "ChiLua/chi_lua.h"

#include "../lbs_linear_boltzmann_solver.h"
#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include <chi_log.h>
extern ChiLog& chi_log;

//###################################################################
/**Initializes the solver.

\param SolverIndex int Handle to the solver.
 \ingroup LuaNPT
 */
int chiLBSInitialize(lua_State *L)
{
  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmann::Solver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    solver = dynamic_cast<LinearBoltzmann::Solver*>(psolver);

    if (not solver)
    {
      chi_log.Log(LOG_ALLERROR) << "chiLBSInitialize: Incorrect solver-type."
                                   " Cannot cast to LinearBoltzmann::Solver\n";
      exit(EXIT_FAILURE);
    }
  }
  catch(const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiLBSInitialize: Invalid handle to solver.\n";
    exit(EXIT_FAILURE);
  }

  solver->Initialize();

  return 0;
}
