#include <ChiLua/chi_lua.h>

#include"../Solver/solver_montecarlon.h"

#include <ChiPhysics/chi_physics.h>

extern ChiPhysics chi_physics_handler;

#include <chi_log.h>

//#############################################################################
/** Executes a MonteCarlon solver.

\param SolverHandle int Handle to the montecarlo solver

\ingroup LuaMonteCarlon
\author Jan*/
int chiMonteCarlonExecute(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiMonteCarlonExecute",1,num_args);

  chi_physics::Solver* solver = nullptr;
  try{
    solver = chi_physics_handler.solver_stack.at(lua_tonumber(L,1));
  }
  catch (std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiMonteCarlonExecute: Invalid solver handle. "
      << lua_tonumber(L,1);
    exit(EXIT_FAILURE);
  }

  if (typeid(*solver) == typeid(chi_montecarlon::Solver))
  {
    auto mcsolver = (chi_montecarlon::Solver*)solver;
    mcsolver->Execute();
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiMonteCarlonExecute: Solver pointed to by solver handle is "
      << " not a MonteCarlo solver.";
    exit(EXIT_FAILURE);
  }
  return 0;
}