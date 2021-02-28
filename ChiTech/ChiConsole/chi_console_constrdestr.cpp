#include "chi_console.h"

//############################################################################# Default constructor
/** Default constructor for the console*/
ChiConsole::ChiConsole() noexcept
{
	//========================================== Initializing console
	this->consoleState = luaL_newstate();
	luaL_openlibs(this->consoleState);
	
	//========================================== Registering functions
	#include"../ChiLua/chi_lua_register.h"

	//========================================== Registering events
	this->InitializeLuaEvent("WM_LBUTTONDOWN");
	this->InitializeLuaEvent("WM_RBUTTONDOWN");
	this->InitializeLuaEvent("WM_MBUTTONDOWN");
	this->InitializeLuaEvent("WM_LBUTTONUP");
	this->InitializeLuaEvent("WM_RBUTTONUP");
	this->InitializeLuaEvent("WM_MBUTTONUP");
	this->InitializeLuaEvent("WM_MOUSEWHEEL");
	this->InitializeLuaEvent("WM_MOUSEMOVE");
	this->InitializeLuaEvent("WM_SHIFTBUTTONDN");
	this->InitializeLuaEvent("WM_SHIFTBUTTONUP");
	this->InitializeLuaEvent("WM_CTLBUTTONDN");
	this->InitializeLuaEvent("WM_CTLBUTTONUP");
	this->InitializeLuaEvent("WM_KEYDN");
	this->InitializeLuaEvent("WM_KEYUP");
	this->InitializeLuaEvent("WM_CHAR");
	this->InitializeLuaEvent("WM_MOUSELEAVE");
	this->InitializeLuaEvent("WM_ENTERDN");
	this->InitializeLuaEvent("WM_ENTERUP");
	this->InitializeLuaEvent("WM_SIZE");
	this->InitializeLuaEvent("WM_CLOSE");

	//========================================== Defining global variables
	char global1[] = "chi_displayStatus=true;";
	luaL_dostring(this->consoleState, global1);

	char global3[] = "chi_programTime=0;";
	luaL_dostring(this->consoleState, global3);

	//========================================== Defining main function
	char initialScript[] =
	  "function main()\n print(\"Empty Main Function\")\nend";
	luaL_dostring(this->consoleState, initialScript);

}


