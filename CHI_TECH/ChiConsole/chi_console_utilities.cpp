//#include <windows.h>
#include <iostream>
#include <cstdio>
#include <sstream>
#include "chi_console.h"



//############################################################################# Initialize event
/** Registers an event for further use.
\author Jan */
void ChiConsole::InitializeLuaEvent(const char* eventTitle)
{
	lua_State* L = this->consoleState;
	CSTEvent defaultEvent;
	lua_newtable(L);

	//===================================================== Pushing status
	lua_pushstring(L, "occured"); lua_pushboolean(L, false);	lua_rawset(L, -3);
	//===================================================== Pushing the integer parameters
	lua_pushstring(L, "iPar0"); lua_pushinteger(L, defaultEvent.iPar[0]);	lua_rawset(L, -3);
	lua_pushstring(L, "iPar1"); lua_pushinteger(L, defaultEvent.iPar[1]);	lua_rawset(L, -3);
	lua_pushstring(L, "iPar2"); lua_pushinteger(L, defaultEvent.iPar[2]);	lua_rawset(L, -3);
	lua_pushstring(L, "iPar3"); lua_pushinteger(L, defaultEvent.iPar[3]);	lua_rawset(L, -3);
	lua_pushstring(L, "iPar4"); lua_pushinteger(L, defaultEvent.iPar[4]);	lua_rawset(L, -3);
	lua_pushstring(L, "iPar5"); lua_pushinteger(L, defaultEvent.iPar[5]);	lua_rawset(L, -3);

	//===================================================== Pushing the boolean parameters
	lua_pushstring(L, "bPar0"); lua_pushboolean(L, defaultEvent.bPar[0]);	lua_rawset(L, -3);
	lua_pushstring(L, "bPar1"); lua_pushboolean(L, defaultEvent.bPar[1]);	lua_rawset(L, -3);
	lua_pushstring(L, "bPar2"); lua_pushboolean(L, defaultEvent.bPar[2]);	lua_rawset(L, -3);
	lua_pushstring(L, "bPar3"); lua_pushboolean(L, defaultEvent.bPar[3]);	lua_rawset(L, -3);
	lua_pushstring(L, "bPar4"); lua_pushboolean(L, defaultEvent.bPar[4]);	lua_rawset(L, -3);
	lua_pushstring(L, "bPar5"); lua_pushboolean(L, defaultEvent.bPar[5]);	lua_rawset(L, -3);

	//===================================================== Pushing the float parameters
	lua_pushstring(L, "fPar0"); lua_pushnumber(L, defaultEvent.fPar[0]);	lua_rawset(L, -3);
	lua_pushstring(L, "fPar1"); lua_pushnumber(L, defaultEvent.fPar[1]);	lua_rawset(L, -3);
	lua_pushstring(L, "fPar2"); lua_pushnumber(L, defaultEvent.fPar[2]);	lua_rawset(L, -3);
	lua_pushstring(L, "fPar3"); lua_pushnumber(L, defaultEvent.fPar[3]);	lua_rawset(L, -3);
	lua_pushstring(L, "fPar4"); lua_pushnumber(L, defaultEvent.fPar[4]);	lua_rawset(L, -3);
	lua_pushstring(L, "fPar5"); lua_pushnumber(L, defaultEvent.fPar[5]);	lua_rawset(L, -3);

	//===================================================== Pushing the strings
	lua_pushstring(L, "sPar0"); lua_pushstring(L, defaultEvent.sPar[0]);	lua_rawset(L, -3);
	lua_pushstring(L, "sPar1"); lua_pushstring(L, defaultEvent.sPar[1]);	lua_rawset(L, -3);
	lua_pushstring(L, "sPar2"); lua_pushstring(L, defaultEvent.sPar[2]);	lua_rawset(L, -3);
	lua_pushstring(L, "sPar3"); lua_pushstring(L, defaultEvent.sPar[3]);	lua_rawset(L, -3);
	lua_pushstring(L, "sPar4"); lua_pushstring(L, defaultEvent.sPar[4]);	lua_rawset(L, -3);
	lua_pushstring(L, "sPar5"); lua_pushstring(L, defaultEvent.sPar[5]);	lua_rawset(L, -3);

	lua_setglobal(L, eventTitle);
	return;
}

//############################################################################# Execute file
/** Executes the given file in the Lua engine.
\author Jan*/
void ChiConsole::ExecuteFile(const char* fileName,int argc, char** argv)
{
	lua_State* L = this->consoleState;
	if (fileName!=NULL)
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
		}
		//std::cout<<lua_tostring(this->consoleState,-1)<<"\n";
	}
}


//############################################################################# Post event to console
/**Posts an event and all of its parameters to the lua-state.

\author CHI Vermaak*/
void ChiConsole::PostEventToConsole(CSTEvent* inputEvent)
{
	lua_State* L = this->consoleState;
	lua_newtable(L);

	//===================================================== Pushing status
	inputEvent->bPar[0]=true;
	lua_pushstring(L, "occured"); lua_pushboolean(L, inputEvent->bPar[0]);	lua_rawset(L, -3);

	//===================================================== Pushing the integer parameters
	lua_pushstring(L, "iPar0"); lua_pushinteger(L, inputEvent->iPar[0]);	lua_rawset(L, -3);
	lua_pushstring(L, "iPar1"); lua_pushinteger(L, inputEvent->iPar[1]);	lua_rawset(L, -3);
	lua_pushstring(L, "iPar2"); lua_pushinteger(L, inputEvent->iPar[2]);	lua_rawset(L, -3);
	lua_pushstring(L, "iPar3"); lua_pushinteger(L, inputEvent->iPar[3]);	lua_rawset(L, -3);
	lua_pushstring(L, "iPar4"); lua_pushinteger(L, inputEvent->iPar[4]);	lua_rawset(L, -3);
	lua_pushstring(L, "iPar5"); lua_pushinteger(L, inputEvent->iPar[5]);	lua_rawset(L, -3);

	//===================================================== Pushing the boolean parameters
	lua_pushstring(L, "bPar0"); lua_pushboolean(L, inputEvent->bPar[0]);	lua_rawset(L, -3);
	lua_pushstring(L, "bPar1"); lua_pushboolean(L, inputEvent->bPar[1]);	lua_rawset(L, -3);
	lua_pushstring(L, "bPar2"); lua_pushboolean(L, inputEvent->bPar[2]);	lua_rawset(L, -3);
	lua_pushstring(L, "bPar3"); lua_pushboolean(L, inputEvent->bPar[3]);	lua_rawset(L, -3);

	//===================================================== Pushing the float parameters
	lua_pushstring(L, "fPar0"); lua_pushnumber(L, inputEvent->fPar[0]);	lua_rawset(L, -3);
	lua_pushstring(L, "fPar1"); lua_pushnumber(L, inputEvent->fPar[1]);	lua_rawset(L, -3);
	lua_pushstring(L, "fPar2"); lua_pushnumber(L, inputEvent->fPar[2]);	lua_rawset(L, -3);
	lua_pushstring(L, "fPar3"); lua_pushnumber(L, inputEvent->fPar[3]);	lua_rawset(L, -3);

	//===================================================== Pushing the strings
	lua_pushstring(L, "sPar0"); lua_pushstring(L, inputEvent->sPar[0]);	lua_rawset(L, -3);
	lua_pushstring(L, "sPar1"); lua_pushstring(L, inputEvent->sPar[1]);	lua_rawset(L, -3);
	lua_pushstring(L, "sPar2"); lua_pushstring(L, inputEvent->sPar[2]);	lua_rawset(L, -3);
	lua_pushstring(L, "sPar3"); lua_pushstring(L, inputEvent->sPar[3]);	lua_rawset(L, -3);



	lua_setglobal(L, inputEvent->eventTitle);
}



//################################################################### Parses MPI
/**Pushes location id and number of processes to lua state.*/
void ChiConsole::PostMPIInfo(int location_id, int number_of_processes)
{
  lua_State* L = this->consoleState;

  lua_pushnumber(L,location_id);
  lua_setglobal(L,"chi_location_id");

  lua_pushnumber(L,number_of_processes);
  lua_setglobal(L,"chi_number_of_processes");
}
