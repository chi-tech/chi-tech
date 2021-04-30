#include "chi_console.h"

#include <iostream>

//############################################################################# Default constructor
/** Default constructor for the console*/
ChiConsole::ChiConsole() noexcept
{
	//========================================== Initializing console
	this->consoleState = luaL_newstate();
	luaL_openlibs(this->consoleState);
	
	//========================================== Registering functions
	#include"../ChiLua/chi_lua_register.h"
}


