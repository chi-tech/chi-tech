#include "../../../CHI_LUA/chi_lua.h"
#include "../chi_thermoalpha.h"
#include<math.h>

#include "../hydro_volume/hydro_volume.h"
#include "../hydro_boundary/hydro_boundary.h"

extern CHI_VECTOR<CHI_THERMOSYSTEM> chithermoSystemStack;

//########################################################## Create boundary condition
/** int chiThermoCreateBC() Creates a hydrodynamic boundary condition.

\params systemHandle int Handle to the system to which the volume belongs.

\return Returns a unique handle for the created boundary condition.

\ingroup LuaThermoalpha
\author Jan*/
int chiThermoCreateBC(lua_State *L)
{
	HYDRO_BOUNDARY* newBoundary = new HYDRO_BOUNDARY;
	
	int systemHandle = lua_tonumber(L,1);
	CHI_THERMOSYSTEM* curSystem = chithermoSystemStack.GetItem(systemHandle);
	
	if (curSystem == NULL) {return 0;}
	
	int index = curSystem->hydroComponentStack.PushItem((HYDRO_COMPONENT*)newBoundary);
	newBoundary->index=index;

	lua_pushnumber(L,index);
	return 1;
}