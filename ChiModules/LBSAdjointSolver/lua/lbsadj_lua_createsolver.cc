#include "ChiLua/chi_lua.h"

#include "LBSAdjointSolver/lbsadj_solver.h"

#include "ChiPhysics/chi_physics.h"

namespace lbs_adjoint
{
namespace lua_utils
{

//###################################################################
/**Creates a lbs_adjoint::AdjointSolver and places it on the physics stack.

\param SolverName string (Optional) Name of the solver.

\return Returns a handle to the solver.*/
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

  auto solver = new lbs_adjoint::AdjointSolver(solver_name);

  auto& physics_handler = ChiPhysics::GetInstance();
  physics_handler.solver_stack.push_back(solver);
  size_t handle = physics_handler.solver_stack.size()-1;

  lua_pushinteger(L,static_cast<lua_Integer>(handle));
  return 1;
}

} //namespace lua_utils
} //namespace lbs_adjoint