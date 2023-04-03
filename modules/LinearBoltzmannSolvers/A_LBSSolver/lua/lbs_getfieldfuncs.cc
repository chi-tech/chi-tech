#include "A_LBSSolver/lbs_solver.h"

#include "chi_runtime.h"

namespace lbs::common_lua_utils
{

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
  const auto& lbs_solver = chi::GetStackItem<lbs::LBSSolver>(chi::solver_stack,
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
          if (pff == lbs_solver.GetFieldFunctions()[ff])
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

}//namespace lbs::common_lua_utils