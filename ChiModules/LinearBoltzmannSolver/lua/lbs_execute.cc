#include "ChiLua/chi_lua.h"
#include "lbs_lua_utils.h"

#include "../lbs_linear_boltzmann_solver.h"

//###################################################################
/**Executes the LBS solver.
\param SolverIndex int Handle to the solver.
 \ingroup LuaLBS
 */
int chiLBSExecute(lua_State *L)
{
  //============================================= Get pointer to solver
  int solver_index = lua_tonumber(L,1);
  auto lbs_solver = LinearBoltzmann::lua_utils::
    GetSolverByHandle(solver_index, __FUNCTION__);

  lbs_solver->Execute();

  return 0;
}
