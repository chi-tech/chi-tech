#include "ChiLua/chi_lua.h"
#include "lbs_lua_utils.h"

//###################################################################
/**Initializes or reinitializes the materials. This normally happens
 * automatically during solver initialization but if the user wants to
 * swap/change XSs during the run then this will allow the material structures
 * to now deal with the new/changed materials.
 *
\param SolverIndex int Handle to the solver maintaining the information.

\ingroup LuaLBS
\author Jan*/
int chiLBSInitializeMaterials(lua_State *L)
{
  const std::string fname = "chiLBSInitializeMaterials";
  const int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);

  auto& lbs_solver = chi::GetStackItem<lbs::SteadyStateSolver>(chi::solver_stack,
                                                               solver_handle,
                                                               fname);

  lbs_solver.InitMaterials();

  return 0;
}