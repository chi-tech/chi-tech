#include "ChiConsole/chi_console.h"

#include "chi_runtime.h"
#include "chi_log.h"

//############################################################################# Execute file
/** Executes the given file in the Lua engine.
\author Jan*/
int chi_objects::ChiConsole::
  ExecuteFile(const std::string& fileName, int argc, char** argv) const
{
	lua_State* L = this->consoleState;
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
		int error = luaL_dofile(this->consoleState,fileName.c_str());

		if (error > 0)
		{
			chi::log.LogAllError() << "LuaError: "
                             << lua_tostring(this->consoleState, -1);
			return EXIT_FAILURE;
		}
	}
	return EXIT_SUCCESS;
}


//################################################################### Parses MPI
/**Pushes location id and number of processes to lua state.*/
void chi_objects::ChiConsole::PostMPIInfo(int location_id, int number_of_processes) const
{
  lua_State* L = this->consoleState;

  lua_pushnumber(L,location_id);
  lua_setglobal(L,"chi_location_id");

  lua_pushnumber(L,number_of_processes);
  lua_setglobal(L,"chi_number_of_processes");
}

//###################################################################
/**Registers a lua function pointing to a c-function.*/
void chi_objects::ChiConsole::
  RegisterFunction(const std::string &string_name, lua_CFunction function)
{
  lua_register(this->consoleState,string_name.c_str(), function);
}
