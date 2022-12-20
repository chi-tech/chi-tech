#include "../lbts_transient_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

#define PropertyArgCntErr(prop_name) \
throw std::logic_error(fname + ": Insufficient amount of arguments, " + \
                       std::to_string(num_args) + ", for property " + \
                       prop_name)

namespace lbs::lbts_lua_utils
{

//###################################################################
/**Sets a property of the LBS-TransientSolver or any derived class.

\param SolverIndex int Handle to the solver for which the set is to be created.
\param PropertyName string String name for a specific property.


##_

##ProperyName
"TIMESTEP"\n
Expects to be followed by a floating point value representing the timestep.\n\n

"TIMESTART"\n
Sets the start time. [Default=0.0]\n\n

"TIMESTOP"\n
Sets the termination time. [Default=0.1]\n\n

"INHIBIT_ADVANCE"\n
Sets a flag that prevents the time from advancing. Can be used to reject a
timestep. [Default=false]\n\n

"VERBOSITY_LEVEL"\n
Sets the verbosity level. Level 0, nothing gets printed. Level 1, basic info.
 [Default=1]\n\n


\author Zachary Hardy*/
int chiLBTSSetProperty(lua_State* L)
{
  const std::string fname = "chiLBTSSetProperty";
  const int num_args = lua_gettop(L);
  if (num_args < 2)
    LuaPostArgAmountError(fname, 2, num_args);

  //============================================= Get the solver
  LuaCheckNilValue(fname, L, 1);
  const int solver_handle = lua_tointeger(L, 1);

  auto& solver = chi::GetStackItem<lbs::TransientSolver>(chi::solver_stack,
                                                         solver_handle,
                                                         fname);

  //============================================= Get the property
  LuaCheckStringValue(fname, L, 2);
  const std::string property = lua_tostring(L, 2);

  if (property == "TIMESTEP")
  {
    if (num_args != 3) PropertyArgCntErr("TIMESTEP");

    LuaCheckNilValue(fname, L, 3);

    const double dt_input = lua_tonumber(L, 3);

    solver.dt = dt_input;

    chi::log.Log() << solver.TextName() << ": dt set to "
                   << std::to_string(dt_input);
  }
  else if (property == "TIMESTART")
  {
    if (num_args != 3) PropertyArgCntErr("TIMESTART");

    LuaCheckNilValue(fname, L, 3);

    const double t_input = lua_tonumber(L, 3);

    solver.t_start = t_input;

    chi::log.Log() << solver.TextName() << ": t_start set to "
                   << std::to_string(t_input);
  }
  else if (property == "TIMESTOP")
  {
    if (num_args != 3) PropertyArgCntErr("TIMESTOP");

    LuaCheckNilValue(fname, L, 3);

    const double t_input = lua_tonumber(L, 3);

    solver.t_final = t_input;

    chi::log.Log() << solver.TextName() << ": t_final set to "
                   << std::to_string(t_input);
  }
  else if (property == "INHIBIT_ADVANCE")
  {
    if (num_args != 3) PropertyArgCntErr("INHIBIT_ADVANCE");

    LuaCheckNilValue(fname, L, 3);

    const bool inhibit_advance = lua_toboolean(L, 3);

    solver.transient_options.inhibit_advance = inhibit_advance;

    chi::log.Log() << solver.TextName() << ": inhibit_advance set to "
                   << std::to_string(inhibit_advance);
  }
  else if (property == "VERBOSITY_LEVEL")
  {
    if (num_args != 3) PropertyArgCntErr("VERBOSITY_LEVEL");

    LuaCheckNilValue(fname, L, 3);

    const int verbosity_level = lua_tointeger(L, 3);

    solver.transient_options.verbosity_level = verbosity_level;

    chi::log.Log() << solver.TextName() << ": verbosity_level set to "
                   << std::to_string(verbosity_level);
  }
  else
    throw std::logic_error(fname + ": unsupported property name \"" +
                           property + "\".");

  return 0;
}

}//namespace lbs::lbts_lua_utils