#include "ChiLua/chi_lua.h"
#include "lbkes_lua_utils.h"

//###################################################################
/**Executes the LBS solver.
\param SolverIndex int Handle to the solver.
 \ingroup LuaLBS
 */
int chiLBKESExecute(lua_State *L)
{
  //============================================= Get pointer to solver
  int solver_index = lua_tonumber(L,1);
  auto lbkes_solver = LinearBoltzmann::k_eigenvalue_lua_utils::
  GetSolverByHandle(solver_index, __FUNCTION__);

  lbkes_solver->Execute();

  return 0;
}
