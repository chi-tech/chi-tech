#include "../../../ChiLua/chi_lua.h"
#include "../chi_thermoalpha.h"
#include<math.h>

#include "../hydro_volume/hydro_volume.h"
#include "../hydro_boundary/hydro_boundary.h"

extern CHI_VECTOR<CHI_THERMOSYSTEM> chithermoSystemStack;

//########################################################## Create volume with coordinates
/** int chiThermoCreateVolumeFromCoordinates() Creates a hydrodynamic volume using a start and end coordinate.

\params systemHandle int Handle to the system to which the volume belongs.
\params point1		 Table Table with fields {x,y,z} containing the start location.
\params point2		 Table Table with fields {x,y,z} containing the end location.

\return Returns a unique handle for the created volume. (int)

\ingroup LuaThermoalpha
\author Jan*/
int chiThermoCreateVolumeFromCoordinates(lua_State *L)
{
	int systemHandle = lua_tonumber(L,1);
	CHI_THERMOSYSTEM* curSystem = chithermoSystemStack.GetItem(systemHandle);
	
	if (curSystem == NULL) {return 0;}
	
	HYDRO_VOLUME* newVolume = new HYDRO_VOLUME;
	
	lua_getfield(L,2,"x");	double x1 = lua_tonumber(L,-1);	lua_pop(L,1);
	lua_getfield(L,2,"y");	double y1 = lua_tonumber(L,-1);	lua_pop(L,1);
	lua_getfield(L,2,"z");	double z1 = lua_tonumber(L,-1);	lua_pop(L,1);
	
	newVolume->location[0] = x1;
	newVolume->location[1] = y1;
	newVolume->location[2] = z1;
	
	lua_getfield(L,3,"x");	double x2 = lua_tonumber(L,-1);	lua_pop(L,1);
	lua_getfield(L,3,"y");	double y2 = lua_tonumber(L,-1);	lua_pop(L,1);
	lua_getfield(L,3,"z");	double z2 = lua_tonumber(L,-1);	lua_pop(L,1);
	
	newVolume->endpoint[0] = x2;
	newVolume->endpoint[1] = y2;
	newVolume->endpoint[2] = z2;
	
	newVolume->L = sqrt((double)((x2-x1)*(x2-x1)+
	                             (y2-y1)*(y2-y1)+
	                             (z2-z1)*(z2-z1)));
	
	newVolume->phi = asin((z2-z1)/newVolume->L);
	newVolume->theta = acos(sqrt((double)((x2-x1)*(x2-x1)+
	                                      (y2-y1)*(y2-y1)))/newVolume->L);
	
	int index = curSystem->hydroComponentStack.PushItem(newVolume);
	newVolume->index=index;
	lua_pushnumber(L,index);
	return 1;
}