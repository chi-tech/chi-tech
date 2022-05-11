#include "ChiLua/chi_lua.h"
#include "physics_lua_utils.h"

#include "chi_log.h"
extern ChiLog& chi_log;

/** \defgroup LuaSolver Solvers
 * \ingroup LuaPhysics*/

//#############################################################################
/** Initializes the solver at the given handle.

\param solver_handle int Handle to the solver.

\ingroup LuaSolver
\author Jan*/
int chiSolverInitialize(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  int solver_handle = lua_tonumber(L, 1);

  auto& solver = chi_physics::lua_utils::GetSolverByHandle(solver_handle,fname);

  solver.Initialize();

  return 0;
}

//#############################################################################
/** Executes the solver at the given handle.

\param solver_handle int Handle to the solver.

\ingroup LuaSolver
\author Jan*/
int chiSolverExecute(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  int solver_handle = lua_tonumber(L, 1);

  auto& solver = chi_physics::lua_utils::GetSolverByHandle(solver_handle,fname);

  solver.Execute();

  return 0;
}

//###################################################################
/** Sets a basic option of a solver.

\param solver_handle int Handle to the reference solver.
\param option_name   string String-name of the option.
\param option_value  varying The value to assign to the option.

\ingroup LuaSolver
\author Jan*/
int chiSolverSetBasicOption(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  LuaCheckIntegerValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);

  const int         solver_handle = lua_tointeger(L,1);
  const std::string option_name   = lua_tostring(L,2);

  auto& solver = chi_physics::lua_utils::GetSolverByHandle(solver_handle, fname);

  try
  {
    auto& option = solver.basic_options[option_name];

    switch (option.Type())
    {
      case chi_data_types::VaryingDataType::VOID:
      case chi_data_types::VaryingDataType::ARBITRARY_BYTES:
        throw std::logic_error("Solver:" + solver.TextName() +
                               " option:" + option_name + " is of invalid type."
                               " This indicates an implementation problem.");
      case chi_data_types::VaryingDataType::STRING:
        LuaCheckStringValue(fname, L, 3);
        option.SetStringValue(lua_tostring(L, 3));
        chi_log.Log() << "Solver:" << solver.TextName()
        << " option:" << option_name
        << " set to " << option.StringValue()
        << ".";
        break;
      case chi_data_types::VaryingDataType::BOOL:
        LuaCheckBoolValue(fname, L, 3);
        option.SetBoolValue(lua_toboolean(L, 3));
        chi_log.Log() << "Solver:" << solver.TextName()
        << " option:" << option_name
        << " set to " << ((option.BoolValue())? "true" : "false")
        << ".";
        break;
      case chi_data_types::VaryingDataType::INTEGER:
        LuaCheckIntegerValue(fname, L, 3);
        option.SetIntegerValue(lua_tointeger(L, 3));
        chi_log.Log() << "Solver:" << solver.TextName()
        << " option:" << option_name
        << " set to " << option.IntegerValue()
        << ".";
        break;
      case chi_data_types::VaryingDataType::FLOAT:
        LuaCheckNumberValue(fname, L, 3);
        option.SetFloatValue(lua_tonumber(L, 3));
        chi_log.Log() << "Solver:" << solver.TextName()
        << " option:" << option_name
        << " set to " << option.FloatValue()
        << ".";
        break;
    }
  }
  catch (const std::out_of_range& oor)
  {
    chi_log.Log(LOG_0ERROR) << fname << ": " << oor.what();
  }

  return 0;
}
