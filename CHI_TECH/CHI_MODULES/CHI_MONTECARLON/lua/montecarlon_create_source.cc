#include"../../../CHI_LUA/chi_lua.h"

#include"../Solver/solver_montecarlon.h"

#include"../../../CHI_PHYSICS/CHI_SOLVER/chi_solver.h"
#include <CHI_MODULES/CHI_MONTECARLON/Source/mc_base_source.h>
#include <CHI_MODULES/CHI_MONTECARLON/Source/BoundarySource/mc_bndry_source.h>


#include <CHI_PHYSICS/chi_physics.h>
#include <chi_log.h>

extern CHI_LOG chi_log;
extern CHI_PHYSICS chi_physics_handler;



//#############################################################################
/** Creates a simple point source at [0 0 0].
 *
\param SolverHandle int Handle to an existing montecarlo solver.
\param SourceType int Source type identifier. See SourceType below.

##_

###PropertyIndex\n
BOUNDARY_SOURCE\n
 Source on a surface boundary. Expects to be followed by the boundary number.\n\n

\return Handle int Handle to the created source.
\ingroup LuaMonteCarlon
\author Jan*/
int chiMonteCarlonCreateSource(lua_State *L)
{
  int num_args = lua_gettop(L);

  int solver_index = lua_tonumber(L,1);
  chi_montecarlon::Solver* solver;

  try{
    solver = (chi_montecarlon::Solver*)chi_physics_handler.solver_stack.at(solver_index);
  }
  catch(std::out_of_range o){
    std::cout << "Invalid solver handle" << std::endl;
    lua_pushinteger(L,-1);
  }

  int source_type = lua_tonumber(L,2);

  //============================================= Boundary source
  if (source_type == MC_BNDRY_SRC)
  {
    if (num_args < 3)
      LuaPostArgAmountError("chiMonteCarlonCreateSource-"
                            "BOUNDARY_SOURCE",
                            4,num_args);

    int ref_boundary = lua_tonumber(L,3);
    if (ref_boundary == 0)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Invalid boundary number supplied in call to "
        << "chiMonteCarlonCreateSource-MC_BNDRY_SRC. Expected a positive number"
           " or a predefined identifier.";
      exit(EXIT_FAILURE);
    }


    chi_montecarlon::BoundarySource* new_source = new
      chi_montecarlon::BoundarySource;

    new_source->ref_bndry = ref_boundary;

    solver->sources.push_back(new_source);
    lua_pushnumber(L,solver->sources.size()-1);

    chi_log.Log(LOG_0) << "MonteCarlo-created boundary source.";
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid boundary type supplied in call to chiMonteCarlonCreateSource";
    exit(EXIT_FAILURE);
  }


  return 1;
}