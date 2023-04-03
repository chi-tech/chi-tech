#include "ChiConsole/chi_console.h"

#include "chi_runtime.h"
#include "chi_log.h"

#define MakeLuaFunctionRegistryEntry(x,y) \
  std::make_pair(x,LuaFunctionRegistryEntry{y,#y})

//############################################################################# Execute file
/** Executes the given file in the Lua engine.
\author Jan*/
int chi_objects::ChiConsole::
  ExecuteFile(const std::string& fileName, int argc, char** argv) const
{
	lua_State* L = this->console_state_;
	if (not fileName.empty())
	{
		if (argc>0)
		{
			lua_newtable(L);
			for (int i=1; i<=argc; i++) {
				lua_pushnumber(L, i);
				lua_pushstring(L, argv[i-1]);
				lua_settable(L, -3);
			}
			lua_setglobal(L,"chiArgs");

		}
		int error = luaL_dofile(this->console_state_, fileName.c_str());

		if (error > 0)
		{
			chi::log.LogAllError() << "LuaError: "
                             << lua_tostring(this->console_state_, -1);
			return EXIT_FAILURE;
		}
	}
	return EXIT_SUCCESS;
}


//###################################################################
/**Pushes location id and number of processes to lua state.*/
void chi_objects::ChiConsole::PostMPIInfo(int location_id, int number_of_processes) const
{
  lua_State* L = this->console_state_;

  lua_pushnumber(L,location_id);
  lua_setglobal(L,"chi_location_id");

  lua_pushnumber(L,number_of_processes);
  lua_setglobal(L,"chi_number_of_processes");
}

//###################################################################
/**Adds a lua_CFunction to the registry. The registry of functions gets
 * parsed into the lua console when `chi::Initialize` is called.*/
char chi_objects::ChiConsole::
  AddFunctionToRegistry(const std::string &raw_name_in_lua,
                        lua_CFunction function_ptr)
{
  //Filter out namespace from the raw name
  std::string name_in_lua = raw_name_in_lua;
  const size_t last_scope = raw_name_in_lua.find_last_of(':');
  if (last_scope != std::string::npos)
    name_in_lua = raw_name_in_lua.substr(last_scope+1, std::string::npos);

  auto& console = GetInstance();

  //Check if the function name is already there
  if (console.lua_function_registry_.count(name_in_lua) > 0)
  {
    const auto& current_entry = console.lua_function_registry_.at(name_in_lua);

    throw std::logic_error(std::string(__PRETTY_FUNCTION__) + ": Attempted "
      "to register lua function \"" + name_in_lua + "\" but the function "
      "is already taken by " + current_entry.function_raw_name);
  }

  console.lua_function_registry_.insert(
    MakeLuaFunctionRegistryEntry(name_in_lua, function_ptr));

  return 0;
}

