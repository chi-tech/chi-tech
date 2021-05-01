#include "chi_log.h"
#include "ChiLua/chi_lua.h"
#include "ChiPhysics/chi_physics.h"

#include "LBSCurvilinear/lbs_curvilinear_solver.h"

extern ChiLog& chi_log;
extern ChiPhysics& chi_physics_handler;


/**Creates a Curvilinear Neutral Particle Transport solver.
 *
 * \return SolverHandle int Handle to the solver created.
 *
 * \code
 * //  cylindrical
 * phys1 = chiLBSCurvilinearCreateSolver(LBSCurvilinear.CYLINDRICAL)
 * //  spherical
 * phys1 = chiLBSCurvilinearCreateSolver(LBSCurvilinear.SPHERICAL)
 * \endcode
 */
int chiLBSCurvilinearCreateSolver(lua_State *L)
{
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);
  const int coord_system = lua_tonumber(L, 1);
  const auto coord_system_type =
    static_cast<chi_math::CoordinateSystemType>(coord_system);

  chi_log.Log(LOG_ALLVERBOSE_1)
    << "Creating Curvilinear Linear Boltzman solver";

  const auto new_solver = new LBSCurvilinear::Solver(coord_system_type);

  chi_physics_handler.solver_stack.push_back(new_solver);
  const int index = chi_physics_handler.solver_stack.size() - 1;
  lua_pushnumber(L,index);

  return 1;
}
