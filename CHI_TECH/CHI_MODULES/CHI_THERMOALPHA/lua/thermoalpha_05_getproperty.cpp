#include "../../../CHI_LUA/chi_lua.h"
#include "../chi_thermoalpha.h"
#include<math.h>

#include "../hydro_volume/hydro_volume.h"
#include "../hydro_boundary/hydro_boundary.h"

extern CHI_VECTOR<CHI_THERMOSYSTEM> chithermoSystemStack;

//######################################################### Set property
/** bool chiThermoGetComponentProperty() Sets the property of a component.
 *
\param sysHndle  Handle to the system being referenced.
\param compHndle Handle to the component being referenced.
\param propCode  Property code.
 
 \return Success=true
\ingroup LuaThermoalpha
\author Jan*/
int chiThermoGetComponentProperty(lua_State *L)
{
	//===================================================== Getting the system
	int systemHandle = lua_tonumber(L,1);
	CHI_THERMOSYSTEM* curSystem = chithermoSystemStack.GetItem(systemHandle);
	
	if (curSystem == NULL) {return false;}
	
	//===================================================== Getting the reference component
	int lcompHandle = lua_tonumber(L,2);
	
	//===================================================== Getting the property code
	int propertyCode = lua_tonumber(L,3);
	
	//===================================================== Processing
	switch (propertyCode)
	{
		case TA_PRESSURE:
		{
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return false;}
			HYDRO_VOLUME* curVolume = (HYDRO_VOLUME*)curComp;
			
			lua_pushnumber(L,curVolume->P);
			return 1;
		}
		case TA_TEMPERATURE:
		{
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return false;}
			HYDRO_VOLUME* curVolume = (HYDRO_VOLUME*)curComp;
			
			lua_pushnumber(L,curVolume->T_f);
			return 1;
		}
		case TA_VELOCITY:
		{
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return false;}
			HYDRO_VOLUME* curVolume = (HYDRO_VOLUME*)curComp;
			
			lua_pushnumber(L,curVolume->v_f);
			return 1;
		}
		case TA_JVELOCITY:
		{
			HYDRO_SJUNCTION* curComp = curSystem->hydroSJunctionStack.GetItem(lcompHandle);
			if (curComp == NULL) {return false;}
			
			lua_pushnumber(L,curComp->v_f);
			return 1;
		}
		case TA_VOIDFRACTION:
		{
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return false;}
			HYDRO_VOLUME* curVolume = (HYDRO_VOLUME*)curComp;
			
			lua_pushnumber(L,curVolume->alpha_g);
			return 1;
		}
		default:
			break;
	}
	return 0;
}