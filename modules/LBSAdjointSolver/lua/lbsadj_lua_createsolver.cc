#include "ChiLua/chi_lua.h"

#include "LBSAdjointSolver/lbsadj_solver.h"

#include "chi_runtime.h"

namespace lbs_adjoint::lua_utils
{

//###################################################################
/**Creates a lbs_adjoint::AdjointSolver and places it on the physics stack.

\param SolverName string (Optional) Name of the solver.

\return Handle Returns a handle to the solver.*/
int chiAdjointSolverCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  std::string solver_name = "AdjointSolver";
  if (num_args == 1)
  {
    LuaCheckNilValue(fname, L, 1);
    LuaCheckStringValue(fname, L, 1);

    solver_name = lua_tostring(L, 1);
  }

  auto solver = std::make_shared<lbs_adjoint::AdjointSolver>(solver_name);

  chi::solver_stack.push_back(solver);
  const size_t handle = chi::solver_stack.size()-1;

  lua_pushinteger(L,static_cast<lua_Integer>(handle));
  return 1;
}

} //namespace lbs_adjoint