#include "CHI_LUA/chi_lua.h"
#include<iostream>
#include "CHI_PHYSICS/chi_physics.h"
#include "CHI_MESH/CHI_REGION/chi_region.h"
#include "CHI_MESH/CHI_MESHHANDLER/chi_meshhandler.h"
#include "CHI_PHYSICS/CHI_FIELDFUNCTION/chi_fieldfunction.h"

#include <chi_log.h>

extern CHI_PHYSICS chi_physics_handler;
extern CHI_LOG     chi_log;


//#############################################################################
/** Adds a field function to a solver.
 *
\param SolverHandle int Handle to the solver.
\param Name char String name for the field function.

\return FFHandle int Global Handle to the field function.

\ingroup LuaSolver
\author Jan*/
int chiSolverAddFieldFunction(lua_State *L)
{
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  int solver_index = lua_tonumber(L,1);
  const char* name = lua_tostring(L,2);

  //======================================================= Getting solver
  chi_physics::Solver* solver;
  try{
    solver = chi_physics_handler.solver_stack.at(solver_index);
  }
  catch(std::out_of_range o)
  {
    chi_log.Log(LOG_ALLERROR)
    << "Invalid solver handle in chiSolverAddFieldFunction";
    exit(EXIT_FAILURE);
  }

  //======================================================= Create Field func
  chi_physics::FieldFunction* new_ff = new chi_physics::FieldFunction;
  new_ff->text_name = std::string(name);
  new_ff->id = chi_physics_handler.fieldfunc_stack.size();

  //======================================================= Add to solver
  solver->field_functions.push_back(new_ff);

  //======================================================= Add to physics
  chi_physics_handler.fieldfunc_stack.push_back(new_ff);
  int index = chi_physics_handler.fieldfunc_stack.size()-1;

  lua_pushnumber(L,index);


  return 1;
}
