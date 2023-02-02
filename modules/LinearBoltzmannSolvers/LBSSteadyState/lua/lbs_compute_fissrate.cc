#include "ChiLua/chi_lua.h"
#include "lbs_lua_utils.h"

//###################################################################
/**Computes and returns the fission rate.
 *
\param SolverIndex int Handle to the solver maintaining the information.
\param OldNewOption string "NEW" or "OLD". For transient solvers the "OLD"
                           option would give the fission rate at the previous
                           timestep. [Default="NEW"]

\return double The fission rate.

\ingroup LuaLBS
\author Jan*/
int chiLBSComputeFissionRate(lua_State *L)
{
  const std::string fname = "chiLBSComputeFissionRate";
  const int num_args = lua_gettop(L);

  if (num_args != 2)
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);

  auto& lbs_solver = chi::GetStackItem<lbs::SteadyStateSolver>(chi::solver_stack,
                                                               solver_handle,
                                                               fname);

  LuaCheckStringValue(fname, L, 2);
  const std::string nature = lua_tostring(L, 2);

  const bool previous = (nature == "OLD");

  const double fission_rate = lbs_solver.ComputeFissionRate(previous);

  lua_pushnumber(L, fission_rate);

  return 1;
}
