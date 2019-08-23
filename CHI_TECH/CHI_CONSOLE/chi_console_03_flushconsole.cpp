#include <iostream>

#include "chi_console.h"
/*


#include "../CHI_GRAPHICS/chi_graphics.h"



extern CHI_GRAPHICS 		chigraphics;
extern CHI_TIMER            chiprogramTimer;
extern CHI_VECTOR<char>     chiconsoleInputBuffer;

*/
#ifdef CHI_USEGRAPHICS
	#include "../CHI_GRAPHICS/chi_graphics.h"
	extern CHI_GRAPHICS 		chigraphics;
#endif
#include "../CHI_LIB/chi_lib.h"
#include "../CHI_PHYSICS/chi_physics.h"
#include "../CHI_TIMER/chi_timer.h"
extern CHI_PHYSICS  		chi_physics_handler;
extern CHI_TIMER        chi_program_timer;

//############################################################################# Flush console
/* This function sends the commands contained in the input buffer to
the lua state from where it is executed.
*/
void CHI_CONSOLE::flushConsole()
{
#ifdef CHI_USEGRAPHICS
	//=============================================================== Uploading events
	int numberOfEvents= this->eventsStack.itemCount;
	for (int k = 0; k < numberOfEvents; k++)
	{
		CST_EVENT* currentEvent = this->eventsStack.PopItem();
		if (currentEvent!= NULL)
		{
			this->PostEventToConsole(currentEvent);
			//try
			//{
			//	delete[] currentEvent;
			//}
//
			//catch (exception& err)
			//{
			//	//MessageBox(NULL,err.what(),"delete Event Exception",MB_OK);
			//}
		}

	}

	//=============================================================== Uploading timestep
	lua_pushnumber(this->consoleState, 16.6667);
	lua_setglobal(this->consoleState, "chi_timestep");

	//=============================================================== Uploading timestep
	lua_pushnumber(this->consoleState, chi_physics_handler.physicsTimeCost);
	lua_setglobal(this->consoleState, "chi_physicsTimeCost");

	//=============================================================== Uploading profileTimes
	lua_pushnumber(this->consoleState, chigraphics.profilerTimes[0]);
	lua_setglobal(this->consoleState, "chi_profileTime0");
	lua_pushnumber(this->consoleState, chigraphics.profilerTimes[1]);
	lua_setglobal(this->consoleState, "chi_profileTime1");
	lua_pushnumber(this->consoleState, chigraphics.profilerTimes[2]);
	lua_setglobal(this->consoleState, "chi_profileTime2");
	lua_pushnumber(this->consoleState, chigraphics.profilerTimes[3]);
	lua_setglobal(this->consoleState, "chi_profileTime3");
	lua_pushnumber(this->consoleState, chigraphics.profilerTimes[4]);
	lua_setglobal(this->consoleState, "chi_profileTime4");
	lua_pushnumber(this->consoleState, chigraphics.profilerTimes[5]);
	lua_setglobal(this->consoleState, "chi_profileTime5");

	//=============================================================== Uploading window size


	//=============================================================== Uploading program time
	lua_pushnumber(this->consoleState, chi_program_timer.GetTime() / 1000.0);
	lua_setglobal(this->consoleState, "chi_programTime");

	//=============================================================== Uploading framerate
	lua_pushnumber(this->consoleState, chigraphics.frameRate);
	lua_setglobal(this->consoleState, "chi_frameRate");

	//=============================================================== Executing the main function
	if (this->runDeveloper)
	{
		luaL_dostring(this->consoleState, "d_main()");
	}
	luaL_dostring(this->consoleState, "main()");



	//=============================================================== Deregister events
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
#endif
	return;
}
