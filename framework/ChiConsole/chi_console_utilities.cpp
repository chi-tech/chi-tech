#include "ChiConsole/chi_console.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_misc_utils.h"

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
  const std::string name_in_lua =
    chi_misc_utils::StringUpToFirstReverse(raw_name_in_lua, "::");

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

//###################################################################
/**Sets/Forms a table structure that mimics the namespace structure of
 * a string. For example the string "sing::sob::nook::Tigger" will be
 * assigned a table structure
 * `sing.sob.nook.Tigger = "sing::sob::nook::Tigger"`.*/
void chi_objects::ChiConsole::SetNamespaceTableStructure(const std::string &key)
{
  auto L = GetInstance().console_state_;
  const auto key_split = chi_misc_utils::StringSplit(key, "::");

  if (key_split.size() == 1)
  {
    lua_pushstring(L, key.c_str());
    lua_setglobal(L, key.c_str());
    return;
  }

  const std::vector<std::string> table_names(key_split.begin(),
                                             key_split.end()-1);

  for (const auto& table_key : table_names)
  {
    // The first entry needs to be in lua's global scope
    // so it looks a little different
    if (table_key == table_names.front())
    {
      lua_getglobal(L, table_names.front().c_str());
      if (not lua_istable(L,-1))
      {
        lua_pop(L,1);
        lua_newtable(L);
        lua_setglobal(L, table_names.front().c_str());
        lua_getglobal(L, table_names.front().c_str());
      }
    }
    else
    {
      lua_getfield(L, -1, table_key.c_str());
      if (not lua_istable(L,-1))
      {
        lua_pop(L,1);
        lua_pushstring(L, table_key.c_str());
        lua_newtable(L);
        lua_settable(L, -3);
        lua_getfield(L, -1, table_key.c_str());
      }
    }

    if (table_key == table_names.back())
    {
      lua_pushstring(L, key_split.back().c_str());
      lua_pushstring(L, key.c_str());
      lua_settable(L,-3);
    }
  }//for table_key in table_keys
  lua_pop(L, lua_gettop(L));
}

