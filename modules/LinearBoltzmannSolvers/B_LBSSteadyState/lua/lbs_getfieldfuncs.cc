#include "ChiLua/chi_lua.h"
#include "lbs_lua_utils.h"

#include "chi_runtime.h"

//###################################################################
/**Obtains a list of field functions from the transport solver.
 *
\param SolverIndex int Handle to the solver for which the list is to be obtained.

\return Pair Table and count. Returns an array of handles and the amount of elements in
        it (indexed from 1).
\ingroup LuaLBS
\author Jan*/
int chiLBSGetFieldFunctionList(lua_State *L)
{
  const std::string fname = "chiLBSGetFieldFunctionList";

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);
  auto& lbs_solver = chi::GetStackItem<lbs::SteadyStateSolver>(chi::solver_stack,
                                                               solver_handle,
                                                               fname);

  //============================================= Push up new table
  lua_newtable(L);
  for (int ff=0; ff<lbs_solver.field_functions.size(); ff++)
  {
    lua_pushnumber(L,ff+1);
    int pff_count = -1;
    bool found = false;
    for (auto& pff : chi::field_function_stack)
    {
      ++pff_count;
      if (pff == lbs_solver.field_functions[ff])
      {
        lua_pushnumber(L,pff_count);
        found = true;
        break;
      }
    }
    if (not found)
      throw std::logic_error(fname + ": SteadyStateSolver field functions not found "
                                     "in global stack.");

    lua_settable(L,-3);
  }

  lua_pushnumber(L,static_cast<lua_Number>(lbs_solver.field_functions.size()));

return 2;
}

//###################################################################
/**Obtains a list of field functions, related only to scalar flux,
from the transport solver.

\param SolverIndex int Handle to the solver for which the list is to be obtained.

\return Pair Table and count. Returns an array of handles and the amount of elements in
        it (indexed from 1).
\ingroup LuaLBS
\author Jan*/
int chiLBSGetScalarFieldFunctionList(lua_State *L)
{
  const std::string fname = "chiLBSGetScalarFieldFunctionList";
  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);
  auto& lbs_solver = chi::GetStackItem<lbs::SteadyStateSolver>(chi::solver_stack,
                                                               solver_handle,
                                                               fname);

  //============================================= Push up new table
  lua_newtable(L);
  int ff=-1;
  int count=0;

  for (int g=0; g<lbs_solver.NumGroups(); g++)
  {
    for (int m=0; m<lbs_solver.NumMoments(); m++)
    {
      ff++;
      if (m==0)
      {
        count++;
        lua_pushnumber(L,count);
        int pff_count = -1;
        bool found = false;
        for (auto& pff : chi::field_function_stack)
        {
          ++pff_count;
          if (pff == lbs_solver.field_functions[ff])
          {
            lua_pushnumber(L,pff_count);
            found = true;
            break;
          }
        }

        if (not found)
          throw std::logic_error(fname + ": SteadyStateSolver field functions not found "
                                         "in global stack.");
        lua_settable(L,-3);
      }
    }
  }

  lua_pushnumber(L,count);
  return 2;
}
