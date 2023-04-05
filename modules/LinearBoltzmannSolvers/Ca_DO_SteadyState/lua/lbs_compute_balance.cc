#include "Ca_DO_SteadyState/lbs_DO_steady_state.h"

#include "chi_runtime.h"
namespace lbs::disc_ord_steady_state_lua_utils
{

//###################################################################
/**Computes balance tables and prints it to the console.
 *
\param SolverIndex int Handle to the solver for which the list is to be obtained.

\ingroup LuaLBS
\author Jan*/
int chiLBSComputeBalance(lua_State *L)
{
  const std::string fname = "chiLBSComputeBalance";
  const int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);

  auto& lbs_solver = chi::GetStackItem<lbs::DiscOrdSteadyStateSolver>(chi::solver_stack,
                                                                      solver_handle,
                                                                      fname);

  lbs_solver.ComputeBalance();

  return 0;
}

}//namespace lbs::disc_ord_steady_state_lua_utils