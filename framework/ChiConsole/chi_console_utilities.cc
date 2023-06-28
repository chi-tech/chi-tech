#include "ChiConsole/chi_console.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_misc_utils.h"

#include "ChiObjectFactory.h"

namespace chi::lua_utils
{
int chiMakeObject(lua_State* L);
}

// #############################################################################
// Execute file
/** Executes the given file in the Lua engine.
\author Jan*/
int chi::ChiConsole::ExecuteFile(const std::string& fileName,
                                 int argc,
                                 char** argv) const
{
  lua_State* L = this->console_state_;
  if (not fileName.empty())
  {
    if (argc > 0)
    {
      lua_newtable(L);
      for (int i = 1; i <= argc; i++)
      {
        lua_pushnumber(L, i);
        lua_pushstring(L, argv[i - 1]);
        lua_settable(L, -3);
      }
      lua_setglobal(L, "chiArgs");
    }
    int error = luaL_dofile(this->console_state_, fileName.c_str());

    if (error > 0)
    {
      Chi::log.LogAllError()
        << "LuaError: " << lua_tostring(this->console_state_, -1);
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}

// ###################################################################
/**Pushes location id and number of processes to lua state.*/
void chi::ChiConsole::PostMPIInfo(int location_id,
                                  int number_of_processes) const
{
  lua_State* L = this->console_state_;

  lua_pushnumber(L, location_id);
  lua_setglobal(L, "chi_location_id");

  lua_pushnumber(L, number_of_processes);
  lua_setglobal(L, "chi_number_of_processes");
}

// ###################################################################
/**Basic addition to registry. Used by the other public methods
 * to registry a text-key to a lua function.*/
void chi::ChiConsole::AddFunctionToRegistry(const std::string& name_in_lua,
                                            lua_CFunction function_ptr)
{
  auto& console = GetInstance();

  // Check if the function name is already there
  if (console.lua_function_registry_.count(name_in_lua) > 0)
  {
    const auto& current_entry = console.lua_function_registry_.at(name_in_lua);

    throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                           ": Attempted "
                           "to register lua function \"" +
                           name_in_lua +
                           "\" but the function "
                           "is already taken by " +
                           current_entry.function_raw_name);
  }

  console.lua_function_registry_.insert(std::make_pair(
    name_in_lua, LuaFunctionRegistryEntry{function_ptr, name_in_lua}));
}

// ###################################################################
/**Adds a lua_CFunction to the registry. The registry of functions gets
 * parsed into the lua console when `chi::Initialize` is called. This
 * particular function will strip the namespace from the the parameter
 * `raw_name_in_lua` and cause the function to be registered in the
 * global namespace of the lua console.*/
char chi::ChiConsole::AddFunctionToRegistryGlobalNamespace(
  const std::string& raw_name_in_lua, lua_CFunction function_ptr)
{
  // Filter out namespace from the raw name
  const std::string name_in_lua =
    chi_misc_utils::StringUpToFirstReverse(raw_name_in_lua, "::");

  AddFunctionToRegistry(name_in_lua, function_ptr);

  return 0;
}

// ###################################################################
/**Adds a lua_CFunction to the registry. The registry of functions gets
 * parsed into the lua console when `chi::Initialize` is called. The full
 * path of the function will be derived from `namespace_name` + "::" +
 * `function_name`.*/
char chi::ChiConsole::AddFunctionToRegistryInNamespaceWithName(
  lua_CFunction function_ptr,
  const std::string& namespace_name,
  const std::string& function_name,
  bool self_callable /*=false*/)
{
  const std::string name_in_lua = namespace_name + "::" + function_name;

  AddFunctionToRegistry(name_in_lua, function_ptr);

  if (self_callable and not namespace_name.empty())
  {
    auto& console = ChiConsole::GetInstance();
    console.class_method_registry_[namespace_name].push_back(function_name);
  }

  return 0;
}

// ###################################################################
chi::InputParameters chi::ChiConsole::DefaultGetInParamsFunc()
{
  return InputParameters();
}

// ###################################################################
/**Wrapper functions operate with input and output parameters, essentially
 * hiding the lua interface.*/
char chi::ChiConsole::AddWrapperToRegistryInNamespaceWithName(
  const std::string& namespace_name,
  const std::string& name_in_lua,
  WrapperGetInParamsFunc syntax_function,
  WrapperCallFunc actual_function,
  bool ignore_null_call_func /*=false*/)
{
  const std::string name = (namespace_name.empty())
                             ? name_in_lua
                             : namespace_name + "::" + name_in_lua;

  auto& console = GetInstance();
  auto& registry = console.function_wrapper_registry_;

  ChiLogicalErrorIf(
    registry.count(name) > 0,
    std::string("Attempted to register lua-function wrapper \"") + name +
      "\" but a wrapper with the same name already exists");

  if (not syntax_function) syntax_function = DefaultGetInParamsFunc;

  if (not ignore_null_call_func)
    ChiLogicalErrorIf(not actual_function, "Problem with get_in_params_func");

  LuaFuncWrapperRegEntry reg_entry;
  reg_entry.get_in_params_func = syntax_function;
  reg_entry.call_func = actual_function;

  registry.insert(std::make_pair(name, reg_entry));

  return 0;
}

// ###################################################################
/**Sets/Forms a lua function in the state using a namespace structure.*/
void chi::ChiConsole::SetLuaFuncNamespaceTableStructure(
  const std::string& full_lua_name, lua_CFunction function_ptr)
{
  auto L = GetInstance().console_state_;
  const auto lua_name_split = chi_misc_utils::StringSplit(full_lua_name, "::");

  if (lua_name_split.size() == 1)
  {
    lua_pushcfunction(L, function_ptr);
    lua_setglobal(L, lua_name_split.back().c_str());
    return;
  }

  const std::vector<std::string> table_names(lua_name_split.begin(),
                                             lua_name_split.end() - 1);

  FleshOutLuaTableStructure(table_names);

  lua_pushstring(L, lua_name_split.back().c_str());
  lua_pushcfunction(L, function_ptr);
  lua_settable(L, -3);

  lua_pop(L, lua_gettop(L));
}

// ###################################################################
/**Sets/Forms a table structure that mimics the namespace structure of
 * a string. For example the string "sing::sob::nook::Tigger" will be
 * assigned a table structure
 * `sing.sob.nook.Tigger = "sing::sob::nook::Tigger"`. Then finally assigns
 * lua call to this table.*/
void chi::ChiConsole::SetLuaFuncWrapperNamespaceTableStructure(
  const std::string& full_lua_name)
{
  auto L = GetInstance().console_state_;

  /**Lambda for making a chunk*/
  auto MakeChunk = [&L, &full_lua_name]()
  {
    std::string chunk_code = "local params = ...; ";
    chunk_code +=
      "return chi_console.LuaWrapperCall(\"" + full_lua_name + "\", ...)";

    luaL_loadstring(L, chunk_code.c_str());
  };

  const auto table_names = chi_misc_utils::StringSplit(full_lua_name, "::");
  std::vector<std::string> namespace_names;
  for (const auto& table_name : table_names)
    if (table_name != table_names.back()) namespace_names.push_back(table_name);

  const auto& function_name = table_names.back();

  if (not namespace_names.empty())
  {
    FleshOutLuaTableStructure(namespace_names);
    lua_pushstring(L, function_name.c_str());
    MakeChunk();
    lua_settable(L, -3);
  }
  else
  {
    MakeChunk();
    lua_setglobal(L, function_name.c_str());
  }

  lua_pop(L, lua_gettop(L));
}

// ###################################################################
/**Sets/Forms a table structure that mimics the namespace structure of
 * a string. For example the string "sing::sob::nook::Tigger" will be
 * assigned a table structure
 * `sing.sob.nook.Tigger = "sing::sob::nook::Tigger"`.*/
void chi::ChiConsole::SetObjectNamespaceTableStructure(
  const std::string& full_lua_name)
{
  auto L = GetInstance().console_state_;

  /**Lambda for registering object type and creation function.*/
  auto RegisterObjectItems = [&L](const std::string& full_name)
  {
    lua_pushstring(L, "type");
    lua_pushstring(L, full_name.c_str());
    lua_settable(L, -3);

    lua_pushstring(L, "Create");
    std::string chunk_code = "local params = ...; ";
    chunk_code += "return chiMakeObjectType(\"" + full_name + "\", ...)";

    luaL_loadstring(L, chunk_code.c_str());
    lua_settable(L, -3);
  };

  const auto table_names = chi_misc_utils::StringSplit(full_lua_name, "::");

  FleshOutLuaTableStructure(table_names);

  RegisterObjectItems(full_lua_name);

  lua_pop(L, lua_gettop(L));
}

// ##################################################################
/**Fleshes out a path in a table tree. For example, given
 * "fee::foo::fah::koo, this routine will make sure that
 * fee.foo.fah.koo is defined as a table tree structure. The routine will
 * create a table structure where one is needed and leave existing ones alone.
 *
 * At the end of the routine the last table in the structure will be on top
 * of the stack.*/
void chi::ChiConsole::FleshOutLuaTableStructure(
  const std::vector<std::string>& table_names)
{
  auto L = GetInstance().console_state_;

  for (const auto& table_name : table_names)
  {
    // The first entry needs to be in lua's global scope,
    // so it looks a little different
    if (table_name == table_names.front())
    {
      lua_getglobal(L, table_name.c_str());
      if (not lua_istable(L, -1))
      {
        lua_pop(L, 1);
        lua_newtable(L);
        lua_setglobal(L, table_name.c_str());
        lua_getglobal(L, table_name.c_str());
      }
    }
    else
    {
      lua_getfield(L, -1, table_name.c_str());
      if (not lua_istable(L, -1))
      {
        lua_pop(L, 1);
        lua_pushstring(L, table_name.c_str());
        lua_newtable(L);
        lua_settable(L, -3);
        lua_getfield(L, -1, table_name.c_str());
      }
    }
  } // for table_key in table_keys
}

// ##################################################################
/**Assumes a table is on top of the stack, then loads the table
 * with chunks that call registered methods of the class with the first
 * argument being the object's handle.*/
void chi::ChiConsole::SetObjectMethodsToTable(const std::string& class_name,
                                              size_t handle)
{
  auto& console = GetInstance();
  auto L = console.console_state_;

  const auto& class_method_registry = console.class_method_registry_;
  if (class_method_registry.count(class_name) == 0) return;

  const auto& method_list = class_method_registry.at(class_name);
  for (const auto& method_name : method_list)
  {
    std::string function_name = class_name;
    function_name += "::";
    function_name += method_name;

    const auto path_names = chi_misc_utils::StringSplit(function_name, "::");
    std::string path_dot_name;
    for (const auto& name : path_names)
    {
      path_dot_name += name;
      if (name != path_names.back()) path_dot_name += ".";
    }

    std::string chunk_code = "local params = ...; ";
    chunk_code +=
      "return " + path_dot_name + "(" + std::to_string(handle) + ", ...)";

    lua_pushstring(L, method_name.c_str());
    luaL_loadstring(L, chunk_code.c_str());
    lua_settable(L, -3);
  }
}

// ##################################################################
/**Makes a formatted output, readible by the documentation scripts,
 * of all the lua wrapper functions.*/
void chi::ChiConsole::DumpRegister() const
{
  Chi::log.Log() << "\n\n";
  for (const auto& [key, entry] : function_wrapper_registry_)
  {
    if (Chi::log.GetVerbosity() == 0)
    {
      Chi::log.Log() << key;
      continue;
    }

    Chi::log.Log() << "LUA_FUNCWRAPPER_BEGIN " << key;

    if (not entry.call_func) Chi::log.Log() << "SYNTAX_BLOCK";

    const auto in_params = entry.get_in_params_func();
    in_params.DumpParameters();

    Chi::log.Log() << "LUA_FUNCWRAPPER_END\n\n";
  }
  Chi::log.Log() << "\n\n";
}

// ##################################################################
/**Given an old status, will update the bindings for only newly registered
 * items.*/
void chi::ChiConsole::UpdateConsoleBindings(
  const chi::RegistryStatuses& old_statuses)
{
  auto ListHasValue =
    [](const std::vector<std::string>& list, const std::string& value)
  { return std::find(list.begin(), list.end(), value) != list.end(); };

  const auto& object_factory = ChiObjectFactory::GetInstance();
  for (const auto& [key, _] : object_factory.Registry())
    if (not ListHasValue(old_statuses.objfactory_keys_, key))
      SetObjectNamespaceTableStructure(key);

  for (const auto& [key, entry] : lua_function_registry_)
    if (not ListHasValue(old_statuses.objfactory_keys_, key))
      SetLuaFuncNamespaceTableStructure(key, entry.function_ptr);

  for (const auto& [key, entry] : function_wrapper_registry_)
    if (not ListHasValue(old_statuses.objfactory_keys_, key))
      if (entry.call_func)
        SetLuaFuncWrapperNamespaceTableStructure(key);
}