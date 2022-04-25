#include <ChiLua/chi_lua.h>

#include "ChiPhysics/chi_physics.h"

#include <chi_log.h>

extern ChiPhysics&  chi_physics_handler;
extern ChiLog     chi_log;

//###################################################################
/**Obtains a named list of the field functions associated with a solver.

\param SolverHandle int A handle to the reference solver.

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
  for (size_t ff=0; ff<solver->field_functions.size(); ff++)
  {
    lua_pushinteger(L,static_cast<lua_Integer>(ff)+1);
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

  lua_pushinteger(L,static_cast<lua_Integer>(solver->field_functions.size()));

  return 2;
}


//###################################################################
/**Gets a field-function handle by name.
\param FFname string Name of the field function.

\return handle If the field-function was found and a handle identified the valid
               handle will be returned (i.e., a natural number >= 0). If the
               field-function by the given name was not found then the function
               will return null.*/
int chiGetFieldFunctionHandleByName(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckStringValue(fname, L ,1);

  const std::string ff_name = lua_tostring(L,1);

  size_t ff_handle_counter = 0;
  std::vector<size_t> handles_that_matched;
  for (const auto& pff : chi_physics_handler.fieldfunc_stack)
  {
    if (pff->text_name == ff_name)
      handles_that_matched.emplace_back(ff_handle_counter);
    ++ff_handle_counter;
  }

  size_t num_handles = handles_that_matched.size();

  if (num_handles == 0)
  {
    chi_log.Log(LOG_0WARNING) << fname << ": No field-functions were found that "
                              << "matched the requested name. A null handle will "
                              << "be returned." << std::endl;

    return 0;
  }

  if (num_handles > 1)
    chi_log.Log(LOG_0WARNING) << fname << ": A total of " << num_handles
                              << " field-functions were found that matched the "
                              << " requested name. Only the first match will be "
                              << " returned.";

  lua_pushinteger(L, static_cast<lua_Integer>(handles_that_matched.front()));
  return 1;
}
