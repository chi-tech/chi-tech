#include "ChiConsole/chi_console.h"

#include <cstdio>

//############################################################################# Execute file
/** Executes the given file in the Lua engine.
\author Jan*/
int chi_objects::ChiConsole::ExecuteFile(const char* fileName, int argc, char** argv) const
{
	lua_State* L = this->consoleState;
	if (fileName != nullptr)
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
		int error = luaL_dofile(this->consoleState,fileName);

		if (error > 0)
		{
			printf("%s\n", lua_tostring(this->consoleState, -1));
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
