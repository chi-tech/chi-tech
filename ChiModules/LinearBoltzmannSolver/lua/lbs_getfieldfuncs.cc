#include "ChiLua/chi_lua.h"
#include "lbs_lua_utils.h"

#include "../lbs_linear_boltzmann_solver.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Obtains a list of field functions from the transport solver.
 *
\param SolverIndex int Handle to the solver for which the list is to be obtained.

\return table,count Returns an array of handles and the amount of elements in
        it (indexed from 1).
\ingroup LuaNPT
\author Jan*/
int chiLBSGetFieldFunctionList(lua_State *L)
{
  //============================================= Get pointer to solver
  int solver_index = lua_tonumber(L,1);
  auto lbs_solver = LinearBoltzmann::lua_utils::
    GetSolverByHandle(solver_index, __FUNCTION__);

  //============================================= Push up new table
  lua_newtable(L);
  for (int ff=0; ff<lbs_solver->field_functions.size(); ff++)
  {
    lua_pushnumber(L,ff+1);
    int pff_count = -1;
    for (auto& pff : chi_physics_handler.fieldfunc_stack)
    {
      ++pff_count;
      if (pff == lbs_solver->field_functions[ff])
      {
        lua_pushnumber(L,pff_count);
        break;
      }
    }

    lua_settable(L,-3);
  }

  lua_pushnumber(L,lbs_solver->field_functions.size());

return 2;
}

//###################################################################
/**Obtains a list of field functions, related only to scalar flux,
from the transport solver.

\param SolverIndex int Handle to the solver for which the list is to be obtained.

\return table,count Returns an array of handles and the amount of elements in
        it (indexed from 1).
\ingroup LuaNPT
\author Jan*/
int chiLBSGetScalarFieldFunctionList(lua_State *L)
{
  //============================================= Get pointer to solver
  int solver_index = lua_tonumber(L,1);
  auto lbs_solver = LinearBoltzmann::lua_utils::GetSolverByHandle(solver_index, "chiLBSGetScalarFieldFunctionList");

  //============================================= Push up new table
  lua_newtable(L);
  int ff=-1;
  int count=0;

  for (int g=0; g<lbs_solver->groups.size(); g++)
  {
    for (int m=0; m<lbs_solver->num_moments; m++)
    {
      ff++;
      if (m==0)
      {
        count++;
        lua_pushnumber(L,count);
        int pff_count = -1;
        for (auto& pff : chi_physics_handler.fieldfunc_stack)
        {
          ++pff_count;
          if (pff == lbs_solver->field_functions[ff])
          {
            lua_pushnumber(L,pff_count);
            break;
          }
        }
        lua_settable(L,-3);
      }
    }
  }

  lua_pushnumber(L,count);
  return 2;
}
