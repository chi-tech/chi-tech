#include"ChiLua/chi_lua.h"

#include"../Solver/solver_montecarlon.h"

#include"ChiPhysics/SolverBase/chi_solver.h"
#include "../Source/mc_base_source.h"
#include "../Source/BoundarySource/mc_bndry_source.h"
#include "../Source/ResidualSource/mc_rmc_source.h"
#include "../Source/ResidualSource/mc_moc_source.h"

#include <ChiPhysics/chi_physics.h>
#include <chi_log.h>

extern ChiLog chi_log;
extern ChiPhysics chi_physics_handler;



//#############################################################################
/** Creates a simple point source at [0 0 0].
 *
\param SolverHandle int Handle to an existing montecarlo solver.
\param SourceType int Source type identifier. See SourceType below.

##_

###PropertyIndex\n
MC_BNDRY_SRC\n
 Source on a surface boundary. Expects to be followed by the boundary number.\n\n

MC_RESID_SRC\n
 Uses a residual source from a field function. Expects to be followed by a
 field function handle. This will sample the residual.\n\n

MC_RESID_SRC_SU\n
 Same as above but will sample the domain uniformly and adjust the weights of
 each individual particle.\n\n

MC_RESID_MOC\n
 Uses a field-function for computing a residual but uses the Method of
 Characteristics to trace the uncollided portions.\n\n

\return Handle int Handle to the created source.
\ingroup LuaMonteCarlon
\author Jan*/
int chiMonteCarlonCreateSource(lua_State *L)
{
  int num_args = lua_gettop(L);

  int solver_index = lua_tonumber(L,1);
  chi_montecarlon::Solver* solver= nullptr;

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
                            3,num_args);

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
  //============================================= Residual source
  else if (source_type == MC_RESID_SRC)
  {
    if (num_args < 3)
      LuaPostArgAmountError("chiMonteCarlonCreateSource-"
                            "MC_RESID_SRC",
                            3,num_args);

    int ff_handle = lua_tonumber(L,3);
    size_t ff_stack_size = chi_physics_handler.fieldfunc_stack.size();


    chi_physics::FieldFunction* ff;
    try {
      ff = chi_physics_handler.fieldfunc_stack.at(ff_handle);
    }

    catch (std::out_of_range o)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Invalid field function handle supplied in call to "
           "chiMonteCarlonCreateSource-MC_RESID_SRC";
      exit(EXIT_FAILURE);
    }

    chi_montecarlon::ResidualSource* new_source = new
      chi_montecarlon::ResidualSource(ff);

    solver->sources.push_back(new_source);
    lua_pushnumber(L,solver->sources.size()-1);

    chi_log.Log(LOG_0) << "MonteCarlo-created residual source.";

  }
  //============================================= Residual source
  //                                              sampled uniformly
  else if (source_type == MC_RESID_SRC_SU)
  {
    if (num_args < 3)
      LuaPostArgAmountError("chiMonteCarlonCreateSource-"
                            "MC_RESID_SRC_SU",
                            3,num_args);

    int ff_handle = lua_tonumber(L,3);
    size_t ff_stack_size = chi_physics_handler.fieldfunc_stack.size();


    chi_physics::FieldFunction* ff;
    try {
      ff = chi_physics_handler.fieldfunc_stack.at(ff_handle);
    }

    catch (std::out_of_range o)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Invalid field function handle supplied in call to "
           "chiMonteCarlonCreateSource-MC_RESID_SRC_SU";
      exit(EXIT_FAILURE);
    }

    chi_montecarlon::ResidualSource* new_source = new
      chi_montecarlon::ResidualSource(ff,true);

    solver->sources.push_back(new_source);
    lua_pushnumber(L,solver->sources.size()-1);

    chi_log.Log(LOG_0) << "MonteCarlo-created residual source.";

  }
  //============================================= Residual source
  //                                              MOC
  else if (source_type == MC_RESID_MOC)
  {
    if (num_args < 3)
      LuaPostArgAmountError("chiMonteCarlonCreateSource-"
                            "MC_RESID_SRC_SU",
                            3,num_args);

    int ff_handle = lua_tonumber(L,3);
    size_t ff_stack_size = chi_physics_handler.fieldfunc_stack.size();


    chi_physics::FieldFunction* ff;
    try {
      ff = chi_physics_handler.fieldfunc_stack.at(ff_handle);
    }

    catch (std::out_of_range o)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Invalid field function handle supplied in call to "
           "chiMonteCarlonCreateSource-MC_RESID_SRC_SU";
      exit(EXIT_FAILURE);
    }

    chi_montecarlon::ResidualMOCSource* new_source = new
      chi_montecarlon::ResidualMOCSource(ff);

    solver->sources.push_back(new_source);
    lua_pushnumber(L,solver->sources.size()-1);

    chi_log.Log(LOG_0) << "MonteCarlo-created residual source.";

  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid boundary type supplied in call to chiMonteCarlonCreateSource";
    exit(EXIT_FAILURE);
  }


  return 1;
}