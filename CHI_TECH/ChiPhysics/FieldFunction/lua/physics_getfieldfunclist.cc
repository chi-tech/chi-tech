#include <ChiLua/chi_lua.h>

#include "ChiPhysics/chi_physics.h"

#include <chi_log.h>

extern ChiPhysics&  chi_physics_handler;
extern ChiLog     chi_log;

//###################################################################
/**Obtains a named list of the field functions associated with a solver.
 *
 * \param SolverHandle int A handle to the reference solver.

 \ingroup LuaSolver */
int chiGetFieldFunctionList(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiGetFieldFunctionList",1,num_args);


  //======================================================= Getting solver
  int solver_index = lua_tonumber(L,1);
  chi_physics::Solver* solver;
  try{
    solver = chi_physics_handler.solver_stack.at(solver_index);
  }
  catch(const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid solver handle in chiGetFieldFunctionList";
    exit(EXIT_FAILURE);
  }

  //============================================= Push up new table
  lua_newtable(L);
  for (int ff=0; ff<solver->field_functions.size(); ff++)
  {
    lua_pushnumber(L,ff+1);
    int pff_count = -1;
    for (auto& pff : chi_physics_handler.fieldfunc_stack)
    {
      ++pff_count;
      if (pff == solver->field_functions[ff])
      {
        lua_pushnumber(L,pff_count);
        break;
      }
    }
    lua_settable(L,-3);
  }

  lua_pushnumber(L,solver->field_functions.size());

  return 2;
}
