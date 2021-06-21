#include "ChiLua/chi_lua.h"

#include "../lbs_linear_boltzmann_solver.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Obtains a list of field functions from the transport solver.
 *
\param SolverIndex int Handle to the solver for which the list is to be obtained.

\ingroup LuaNPT
\author Jan*/
int chiLBSComputeBalance(lua_State *L)
{
  int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmann::Solver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    solver = dynamic_cast<LinearBoltzmann::Solver*>(psolver);

    if (not psolver)
    {
      chi_log.Log(LOG_ALLERROR) << "chiLBSComputeBalance: Incorrect solver-type."
                                   " Cannot cast to LinearBoltzmann::Solver\n";
      exit(EXIT_FAILURE);
    }
  }
  catch(const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      <<"Invalid handle to solver"
        "in chiLBSComputeBalance\n";
    exit(EXIT_FAILURE);
  }

  solver->ComputeBalance();

  return 0;
}