#include <iostream>

#include "chi_console.h"


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

	return;
}
