#include "ChiLua/chi_lua.h"
#include "lbs_lua_utils.h"

#include "../lbs_linear_boltzmann_solver.h"

//###################################################################
/**Obtains a list of field functions from the transport solver.
 *
\param SolverIndex int Handle to the solver for which the list is to be obtained.

\ingroup LuaLBS
\author Jan*/
int chiLBSComputeBalance(lua_State *L)
{
  int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  //============================================= Get pointer to solver
  int solver_handle = lua_tonumber(L, 1);
  auto& lbs_solver = lbs::lua_utils::
    GetSolverByHandle(solver_handle, __FUNCTION__);

  lbs_solver.ComputeBalance();

  return 0;
}
