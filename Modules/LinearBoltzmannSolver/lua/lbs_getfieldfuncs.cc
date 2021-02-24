#include "ChiLua/chi_lua.h"

#include "../lbs_linear_boltzmann_solver.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include <chi_log.h>
extern ChiLog& chi_log;

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
  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmann::Solver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    solver = dynamic_cast<LinearBoltzmann::Solver*>(psolver);

    if (not psolver)
    {
      chi_log.Log(LOG_ALLERROR) << "chiLBSGetFieldFunctionList: Incorrect solver-type."
                                   " Cannot cast to LinearBoltzmann::Solver\n";
      exit(EXIT_FAILURE);
    }
  }
  catch(const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      <<"Invalid handle to solver"
        "in chiLBSGetFieldFunctionList\n";
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
  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzmann::Solver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    solver = dynamic_cast<LinearBoltzmann::Solver*>(psolver);

    if (not psolver)
    {
      chi_log.Log(LOG_ALLERROR) << "chiLBSGetScalarFieldFunctionList: Incorrect solver-type."
                                   " Cannot cast to LinearBoltzmann::Solver\n";
      exit(EXIT_FAILURE);
    }
  }
  catch(const std::out_of_range& o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiLBSGetScalarFieldFunctionList\n");
    exit(EXIT_FAILURE);
  }

  //============================================= Push up new table
  lua_newtable(L);
  int ff=-1;
  int count=0;

  for (int g=0; g<solver->groups.size(); g++)
  {
    for (int m=0; m<solver->num_moments; m++)
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
          if (pff == solver->field_functions[ff])
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
