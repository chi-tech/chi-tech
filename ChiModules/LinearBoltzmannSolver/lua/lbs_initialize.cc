#include "ChiLua/chi_lua.h"
#include "lbs_lua_utils.h"

#include "../lbs_linear_boltzmann_solver.h"

//###################################################################
/**Initializes the solver.

\param SolverIndex int Handle to the solver.
 \ingroup LuaLBS
 */
int chiLBSInitialize(lua_State *L)
{
  //============================================= Get pointer to solver
  int solver_index = lua_tonumber(L,1);
  auto lbs_solver = LinearBoltzmann::lua_utils::
    GetSolverByHandle(solver_index, __FUNCTION__);

  lbs_solver->Initialize();

  return 0;
}
