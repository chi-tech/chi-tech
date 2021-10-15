#include "ChiLua/chi_lua.h"
#include "lbkes_lua_utils.h"

using namespace LinearBoltzmann;

//###################################################################
/**Initialize the solver.*/
int chiLBKESInitialize(lua_State* L)
{
  //============================================= Get pointer to solver
  int solver_index = lua_tonumber(L,1);
  auto lbkes_solver = LinearBoltzmann::k_eigenvalue_lua_utils::
    GetSolverByHandle(solver_index, __FUNCTION__);

  lbkes_solver->Initialize();

  return 0;
}
