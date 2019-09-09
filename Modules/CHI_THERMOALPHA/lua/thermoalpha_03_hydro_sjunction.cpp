#include "ChiLua/chi_lua.h"
#include "../chi_thermoalpha.h"
#include<math.h>

#include "../hydro_volume/hydro_volume.h"
#include "../hydro_boundary/hydro_boundary.h"

extern CHI_VECTOR<CHI_THERMOSYSTEM> chithermoSystemStack;

//########################################################## Create Single Junction
/** int chiThermoCreateSJunction() Creates a hydrodynamic single junction.

\params systemHandle Handle to the system to which the volume belongs.

\return Returns a unique handle for the created single junction.

\ingroup LuaThermoalpha
\author Jan*/
int chiThermoCreateSJunction(lua_State *L)
{
	HYDRO_SJUNCTION* newSJunction = new HYDRO_SJUNCTION;
	
	int systemHandle = lua_tonumber(L,1);
	CHI_THERMOSYSTEM* curSystem = chithermoSystemStack.GetItem(systemHandle);
	
	if (curSystem == NULL) {return 0;}
	
	int index = curSystem->hydroSJunctionStack.PushItem(newSJunction);
	
	lua_pushnumber(L,index);
	return 1;
}