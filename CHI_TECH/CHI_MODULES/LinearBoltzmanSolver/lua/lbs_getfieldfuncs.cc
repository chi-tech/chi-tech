#include "../../../CHI_LUA/chi_lua.h"

#include "CHI_MODULES/LinearBoltzmanSolver/lbs_linear_boltzman_solver.h"

#include "../../../CHI_PHYSICS/chi_physics.h"
extern CHI_PHYSICS chi_physics_handler;

#include <chi_log.h>
extern CHI_LOG chi_log;

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
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
      <<"ERROR: Incorrect solver-type"
        "in chiLBSSetProperty\n";
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
      <<"ERROR: Invalid handle to solver"
        "in chiLBSSetProperty\n";
    exit(EXIT_FAILURE);
  }

  //============================================= Push up new table
  lua_newtable(L);
  for (int ff=0; ff<solver->field_functions.size(); ff++)
  {
    lua_pushnumber(L,ff+1);
    lua_pushnumber(L,solver->field_functions[ff]->id);
    lua_settable(L,-3);

//    chi_log.Log(LOG_0)
//    << solver->field_functions[ff]->text_name;
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
  LinearBoltzmanSolver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);

    if (typeid(*psolver) == typeid(LinearBoltzmanSolver))
    {
      solver = (LinearBoltzmanSolver*)(psolver);
    }
    else
    {
      fprintf(stderr,"ERROR: Incorrect solver-type"
                     "in chiLBSSetProperty\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(std::out_of_range o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiLBSSetProperty\n");
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
        lua_pushnumber(L,solver->field_functions[ff]->id);
        lua_settable(L,-3);
      }
    }
  }


  lua_pushnumber(L,count);
  return 2;
}