#include "ChiLua/chi_lua.h"

#include "../k_eigenvalue_solver.h"
#include "lua_utils.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics& chi_physics_handler;

#include "chi_log.h"
extern ChiLog& chi_log;


//###################################################################
/**Initializes the k-eigenvalue solver.*/
int chiKEigenvalueLBSInitialize(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  LuaCheckIntegerValue(__FUNCTION__, L, 1);

  int solver_index = lua_tonumber(L, 1);

  //============================================= Get pointer to solver
  auto solver = LinearBoltzmann::lua_utils::
  GetKEigenvalueSolverByHandle(solver_index, __FUNCTION__);

  solver->InitializeKSolver();

  return 1;
}
