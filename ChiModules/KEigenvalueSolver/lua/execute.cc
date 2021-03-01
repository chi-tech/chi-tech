#include "ChiLua/chi_lua.h"

#include "../k_eigenvalue_solver.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include "chi_log.h"
extern ChiLog& chi_log;

using namespace LinearBoltzmann;

//###################################################################
/**Executes the k-eigenvalue solver.*/
int chiKEigenvalueLBSExecute(lua_State *L)
{
  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  KEigenvalue::Solver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    solver = dynamic_cast<KEigenvalue::Solver*>(psolver);

    if (not solver)
    {
      chi_log.Log(LOG_ALLERROR) << "chiLBSKEigenvalueExecute: Incorrect solver-type."
                                   " Cannot cast to KEigenvalue::Solver\n";
      exit(EXIT_FAILURE);
    }
  }
  catch(const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR) << "chiLBSKEigenvalueExecute: Invalid handle to solver\n";
    exit(EXIT_FAILURE);
  }

  solver->ExecuteKSolver();

  return 1;
}

