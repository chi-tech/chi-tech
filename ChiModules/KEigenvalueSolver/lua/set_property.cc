#include "ChiLua/chi_lua.h"

#include "KEigenvalueSolver/k_eigenvalue_solver.h"

#include "ChiPhysics/chi_physics.h"

extern ChiPhysics& chi_physics_handler;

#include <chi_log.h>

extern ChiLog& chi_log;

using namespace LinearBoltzmann;

//############################################################
/**Set boolean for using precursors.*/
int chiLBSSetUsePrecursors(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args != 2)
    LuaPostArgAmountError(__FUNCTION__, 2, num_args);
  LuaCheckNilValue(__FUNCTION__, L, 2);

  int solver_index = lua_tonumber(L, 1);
  bool use_precursors = lua_toboolean(L, 2);

  // ----- Get pointer to solver
  LinearBoltzmann::Solver* solver;
  try
  {

    solver = dynamic_cast<LinearBoltzmann::Solver*>(chi_physics_handler.solver_stack.at(solver_index));

    if (not solver)
    {
      chi_log.Log(LOG_ALLERROR)
          << __FUNCTION__ << ": Incorrect solver-type."
                             " Cannot cast to KEigenvalue::Solver\n";
      exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
        << __FUNCTION__ << ": Invalid handle to solver";
    exit(EXIT_FAILURE);
  }
  solver->options.use_precursors = use_precursors;

  chi_log.Log(LOG_0)
      << "Precursors set to  " << use_precursors;

  return 1;
}

//############################################################
/**Set maximum number of k-eigenvalue iterations.*/
int chiLBSSetMaxKIterations(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args != 2)
    LuaPostArgAmountError(__FUNCTION__, 2, num_args);
  LuaCheckNilValue(__FUNCTION__, L, 2);

  int solver_index = lua_tonumber(L, 1);
  int num_iter = lua_tointeger(L, 2);

  // ----- Get pointer to solver
  chi_physics::Solver* psolver;
  KEigenvalue::Solver* solver;
  try
  {

    solver = dynamic_cast<KEigenvalue::Solver*>(psolver);

    if (not solver)
    {
      chi_log.Log(LOG_ALLERROR)
          << __FUNCTION__ << ": Incorrect solver-type."
                             " Cannot cast to KEigenvalue::Solver\n";
      exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
        << __FUNCTION__ << ": Invalid handle to solver";
    exit(EXIT_FAILURE);
  }

  if (num_iter < 0)
  {
    chi_log.Log(LOG_ALLERROR)
        << "Invalid number of iterations "
        << "in call to chiLBSSetMaxKIterations. "
        << "Must be >= 0.";
    exit(EXIT_FAILURE);
  }
  solver->max_iterations = num_iter;

  chi_log.Log(LOG_0)
      << "Eigenvalue max # iterations set to " << num_iter;

  return 0;
}

//############################################################
/**Set the k-eigenvalue solver tolerance.*/
int chiLBSSetKTolerance(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args != 2)
    LuaPostArgAmountError(__FUNCTION__, 2, num_args);
  LuaCheckNilValue(__FUNCTION__, L, 2);

  int solver_index = lua_tonumber(L, 1);
  double tol = lua_tonumber(L, 2);

  // ----- Get pointer to solver
  chi_physics::Solver* psolver;
  KEigenvalue::Solver* solver;
  try
  {

    solver = dynamic_cast<KEigenvalue::Solver*>(psolver);

    if (not solver)
    {
      chi_log.Log(LOG_ALLERROR)
          << __FUNCTION__ << ": Incorrect solver-type."
                             " Cannot cast to KEigenvalue::Solver\n";
      exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
        << __FUNCTION__ << ": Invalid handle to solver";
    exit(EXIT_FAILURE);
  }

  if (tol < 0)
  {
    chi_log.Log(LOG_ALLERROR)
        << "Invalid tolerance specified. Must be >= 0.0";
    exit(EXIT_FAILURE);
  }
  solver->tolerance = tol;

  char buff[100];
  sprintf(buff, "%.4e", tol);

  chi_log.Log(LOG_0)
      << "k-eigenvalue tolerance set to " << buff;

  return 0;
}
