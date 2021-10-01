#include "../lbkes_k_eigenvalue_solver.h"

#include <chi_lua.h>

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics& chi_physics_handler;

#include <chi_log.h>
extern ChiLog& chi_log;

using namespace LinearBoltzmann;

//###################################################################
/**Executes the solver.*/
int chiLBKESExecute(lua_State* L)
{
  int solver_index = lua_tonumber(L, 1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  KEigenvalueSolver* solver;
  try
  {
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    solver = dynamic_cast<KEigenvalueSolver*>(psolver);

    if (not solver)
    {
      chi_log.Log(LOG_ALLERROR)
          << __FUNCTION__ << ": Incorrect solver-type. "
          << "Cannot cast to LinearBoltzmann::KEigenvalueSolver.";
      exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
        << __FUNCTION__ << ": Invalid handle to solver.";
    exit(EXIT_FAILURE);
  }

  //============================================= Execute
  solver->Execute();

  return 1;
}
