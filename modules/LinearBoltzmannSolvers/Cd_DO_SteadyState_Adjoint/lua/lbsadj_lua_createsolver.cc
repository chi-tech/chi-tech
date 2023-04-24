#include "ChiLua/chi_lua.h"

#include "Cd_DO_SteadyState_Adjoint/lbsadj_solver.h"

#include "chi_runtime.h"

namespace lbs::adjoint_lua_utils
{

//###################################################################
/**Creates a lbs::SteadyStateAdjointSolver and places it on the physics stack.

\param SolverName string (Optional) Name of the solver.

\return Handle Returns a handle to the solver.*/
int chiAdjointSolverCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  std::string solver_name = "SteadyStateAdjointSolver";
  if (num_args == 1)
  {
    LuaCheckNilValue(fname, L, 1);
    LuaCheckStringValue(fname, L, 1);

    solver_name = lua_tostring(L, 1);
  }

  auto solver = std::make_shared<lbs::DiscOrdSteadyStateAdjointSolver>(solver_name);

  chi::object_stack.push_back(solver);
  const size_t handle = chi::object_stack.size()-1;

  lua_pushinteger(L,static_cast<lua_Integer>(handle));
  return 1;
}

} //namespace lbs