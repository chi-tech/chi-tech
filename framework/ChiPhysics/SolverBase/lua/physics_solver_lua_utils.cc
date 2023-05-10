#include "physics_solver_lua_utils.h"

#include "ChiPhysics/SolverBase/chi_solver.h"
#include "ChiPhysics/FieldFunction/fieldfunction_gridbased.h"

#include "ChiObject/object_maker.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiConsole/chi_console.h"

/** \defgroup LuaSolver Solvers
 * \ingroup LuaPhysics*/

namespace chi_physics::lua_utils
{
RegisterLuaFunctionAsIs(chiSolverCreate);

RegisterLuaFunctionAsIs(chiSolverInitialize);
RegisterLuaFunctionAsIs(chiSolverExecute);
RegisterLuaFunctionAsIs(chiSolverStep);
RegisterLuaFunctionAsIs(chiSolverAdvance);
RegisterLuaFunctionAsIs(chiSolverSetBasicOption);
RegisterLuaFunctionAsIs(chiSolverGetName);
RegisterLuaFunctionAsIs(chiSolverGetFieldFunctionList);

// #############################################################################
/**Generic lua routine for the creation of solvers.
 * \param params ParameterBlock. A single block with at least one field
 *                   \"type\", which contains a registered solver type.
 * ## _
 *
 * Example:
\code
chiSolverCreate({type=cfem_diffusion.Solver})
\endcode*/
int chiSolverCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckTableValue(fname, L, 1);

  const auto params = chi_lua::TableParserAsParameterBlock::ParseTable(L, 1);

  const auto& object_maker = ChiObjectMaker::GetInstance();
  const size_t handle = object_maker.MakeRegisteredObject(params);

  lua_pushinteger(L, static_cast<lua_Integer>(handle));
  return 1;
}

// #############################################################################
/** Initializes the solver at the given handle.

\param solver_handle int Handle to the solver.

\ingroup LuaSolver
\author Jan*/
int chiSolverInitialize(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const int solver_handle = lua_tonumber(L, 1);

  auto& solver = chi::GetStackItem<chi_physics::Solver>(
    chi::object_stack, solver_handle, fname);

  solver.Initialize();

  return 0;
}

// #############################################################################
/** Executes the solver at the given handle.

\param solver_handle int Handle to the solver.

\ingroup LuaSolver
\author Jan*/
int chiSolverExecute(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const int solver_handle = lua_tonumber(L, 1);

  auto& solver = chi::GetStackItem<chi_physics::Solver>(
    chi::object_stack, solver_handle, fname);

  solver.Execute();

  return 0;
}

// #############################################################################
/** Performs a single timestep for the solver at the given handle.

\param solver_handle int Handle to the solver.

\ingroup LuaSolver
\author Jan*/
int chiSolverStep(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const int solver_handle = lua_tonumber(L, 1);

  auto& solver = chi::GetStackItem<chi_physics::Solver>(
    chi::object_stack, solver_handle, fname);

  solver.Step();

  return 0;
}

// #############################################################################
/** Advances the time values of the solver at the given handle.

\param solver_handle int Handle to the solver.

\ingroup LuaSolver
\author Jan*/
int chiSolverAdvance(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const int solver_handle = lua_tonumber(L, 1);

  auto& solver = chi::GetStackItem<chi_physics::Solver>(
    chi::object_stack, solver_handle, fname);

  solver.Advance();

  return 0;
}

// ###################################################################
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
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  LuaCheckIntegerValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);

  const int solver_handle = lua_tointeger(L, 1);
  const std::string option_name = lua_tostring(L, 2);

  auto& solver = chi::GetStackItem<chi_physics::Solver>(
    chi::object_stack, solver_handle, fname);

  try
  {
    auto& option = solver.GetBasicOptions()[option_name];

    switch (option.Type())
    {
      case chi_data_types::VaryingDataType::VOID:
      case chi_data_types::VaryingDataType::ARBITRARY_BYTES:
        throw std::logic_error("Solver:" + solver.TextName() +
                               " option:" + option_name +
                               " is of invalid type."
                               " This indicates an implementation problem.");
      case chi_data_types::VaryingDataType::STRING:
        LuaCheckStringValue(fname, L, 3);
        option.SetStringValue(lua_tostring(L, 3));
        chi::log.Log() << "Solver:" << solver.TextName()
                       << " option:" << option_name << " set to "
                       << option.StringValue() << ".";
        break;
      case chi_data_types::VaryingDataType::BOOL:
        LuaCheckBoolValue(fname, L, 3);
        option.SetBoolValue(lua_toboolean(L, 3));
        chi::log.Log() << "Solver:" << solver.TextName()
                       << " option:" << option_name << " set to "
                       << ((option.BoolValue()) ? "true" : "false") << ".";
        break;
      case chi_data_types::VaryingDataType::INTEGER:
        LuaCheckIntegerValue(fname, L, 3);
        option.SetIntegerValue(lua_tointeger(L, 3));
        chi::log.Log() << "Solver:" << solver.TextName()
                       << " option:" << option_name << " set to "
                       << option.IntegerValue() << ".";
        break;
      case chi_data_types::VaryingDataType::FLOAT:
        LuaCheckNumberValue(fname, L, 3);
        option.SetFloatValue(lua_tonumber(L, 3));
        chi::log.Log() << "Solver:" << solver.TextName()
                       << " option:" << option_name << " set to "
                       << option.FloatValue() << ".";
        break;
    }
  }
  catch (const std::out_of_range& oor)
  {
    chi::log.Log0Error() << fname << ": " << oor.what();
    throw oor;
  }

  return 0;
}

// #############################################################################
/** Returns the text name of the solver.

\param solver_handle int Handle to the solver.

\ingroup LuaSolver
\author Jan*/
int chiSolverGetName(lua_State* L)
{
  const std::string fname = "chiSolverGetName";
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);
  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const int solver_handle = lua_tonumber(L, 1);

  const auto& solver = chi::GetStackItem<chi_physics::Solver>(
    chi::object_stack, solver_handle, fname);

  lua_pushstring(L, solver.TextName().c_str());

  return 1;
}

// ###################################################################
/**Obtains a named list of the field functions associated with a solver.

\param SolverHandle int A handle to the reference solver.

\ingroup LuaSolver */
int chiSolverGetFieldFunctionList(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiGetFieldFunctionList", 1, num_args);

  //======================================================= Getting solver
  const int solver_handle = lua_tonumber(L, 1);

  const auto& solver = chi::GetStackItem<chi_physics::Solver>(
    chi::object_stack, solver_handle, fname);

  //============================================= Push up new table
  lua_newtable(L);
  for (size_t ff = 0; ff < solver.GetFieldFunctions().size(); ff++)
  {
    lua_pushinteger(L, static_cast<lua_Integer>(ff) + 1);
    int pff_count = -1;
    bool found = false;
    for (auto& pff : chi::field_function_stack) // pff pointer to field func
    {
      //const auto& compare_ff_ptr = solver.GetFieldFunctions()[ff];
      //const auto compare_base_ptr =
      //  std::static_pointer_cast<chi_physics::FieldFunction>(compare_ff_ptr);
      ++pff_count;
      if (pff == solver.GetFieldFunctions()[ff])
      {
        lua_pushnumber(L, pff_count);
        found = true;
        break;
      }
    }

    if (not found)
      throw std::logic_error(fname + ": The solver specified has no "
                                     "field functions that match the global"
                                     " stack.");
    lua_settable(L, -3);
  }

  lua_pushinteger(L,
                  static_cast<lua_Integer>(solver.GetFieldFunctions().size()));

  return 2;
}
//}//namespace temp
} // namespace chi_physics::lua_utils