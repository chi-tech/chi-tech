#include "chi_physics.h"
#include"../ChiTimer/chi_timer.h"

#include "chi_runtime.h"

#include"../ChiConsole/chi_console.h"
extern ChiConsole&  chi_console;

#include <iostream>
#include <unistd.h>

//############################################################################# Run physics loop
/** Timed loop executing all physics events.*/
void ChiPhysics::RunPhysicsLoop()
{
	ChiTimer physicsTiming;
	ChiTimer profilingTimer;
	double    physicsTime=0.0;
	double    physicsOldTime=0.0;
	int       cycleCollected=0;
	double    timeCollected=0;
	double    profileTime=0.0;


	physicsTiming.Reset();
	profilingTimer.Reset();
	std::cout<<"Physics loop started.\n";

	while (!chi::run_time::termination_posted)
	{
		physicsTime=physicsTiming.GetTime();

		if (physicsTime>=this->physicsTimestep+physicsOldTime)
		{
			profilingTimer.Reset();
			physicsOldTime=physicsTime;

			//============================================= Flushes the console
      chi_console.FlushConsole();
//			for (int k=0;k<chiconsoleInputBuffer.itemCount;k++)
//			{
//         int error = luaL_dostring(chi_console.consoleState, chiconsoleInputBuffer.GetItem(k));
//         if (error > 0)
//         {
//             printf("%s\n", lua_tostring(chi_console.consoleState, -1));
//         }
//			}
//			chiconsoleInputBuffer.ClearVector();

			cycleCollected++;
			timeCollected = timeCollected+profilingTimer.GetTime();
			profilingTimer.Reset();

			if (cycleCollected>=10)
			{
				timeCollected=0;
				cycleCollected=0;
			}
		}

		profileTime = profilingTimer.GetTime();
		usleep(    std::max(0.0,(double)(16666-profileTime*1000)));
	}

	std::cout<<"Physics loop stopped successfully.\n";
}

