#include "../../../CHI_LUA/chi_lua.h"
#include "../chi_thermoalpha.h"
#include<math.h>

#include "../hydro_volume/hydro_volume.h"
#include "../hydro_boundary/hydro_boundary.h"

extern CHI_VECTOR<CHI_THERMOSYSTEM> chithermoSystemStack;

//######################################################### Set property
/** bool chiThermoSetComponentProperty() Sets the property of a component.
 *
\param sysHndle  Handle to the system being referenced.
\param compHndle Handle to the component being referenced.
\param propCode  Property code.
 
 \return Success=true
\ingroup LuaThermoalpha
\author Jan*/
int chiThermoSetComponentProperty(lua_State *L)
{
	//===================================================== Getting the system
	int systemHandle = lua_tonumber(L,1);
	CHI_THERMOSYSTEM* curSystem = chithermoSystemStack.GetItem(systemHandle);
	
	if (curSystem == NULL)
	{
		printf("ERROR: Could not find the specified system!\n");
		return false;
	}
	
	//===================================================== Getting the reference component
	int lcompHandle = lua_tonumber(L,2);
	
	//===================================================== Getting the property code
	int propertyCode = lua_tonumber(L,3);
	
	//===================================================== Processing
	switch (propertyCode)
	{
		case TA_PRESSURE:
		{
			double value = lua_tonumber(L,4);
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return 0;}
			HYDRO_VOLUME* curVolume = (HYDRO_VOLUME*)curComp;
			
			curVolume->P = value;
			break;
		}
		case TA_BCPRESSURE:
		{
			double value = lua_tonumber(L,4);
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return 0;}
			HYDRO_BOUNDARY* curVolume = (HYDRO_BOUNDARY*)curComp;
			
			curVolume->P = value;
			break;
		}
		case TA_TEMPERATURE:
		{
			double value = lua_tonumber(L,4);
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return 0;}
			HYDRO_VOLUME* curVolume = (HYDRO_VOLUME*)curComp;
			
			curVolume->T_f = value;
			curVolume->T_g = value;
			break;
		}
		case TA_BCTEMPERATURE:
		{
			double value = lua_tonumber(L,4);
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return 0;}
			HYDRO_BOUNDARY* curVolume = (HYDRO_BOUNDARY*)curComp;
			
			curVolume->T_f = value;
			curVolume->T_g = value;
			break;
		}
		case TA_VELOCITY:
		{
			double value = lua_tonumber(L,4);
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return 0;}
			HYDRO_VOLUME* curVolume = (HYDRO_VOLUME*)curComp;
			
			curVolume->v_f = value;
			curVolume->v_g = value;
			break;
		}
		case TA_JVELOCITY:
		{
			double value = lua_tonumber(L,4);
			HYDRO_SJUNCTION* curComp = curSystem->hydroSJunctionStack.GetItem(lcompHandle);
			
			curComp->v_f = value;
			curComp->v_g = 0.0;
			break;
		}
		case TA_VOIDFRACTION:
		{
			double value = lua_tonumber(L,4);
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return 0;}
			HYDRO_VOLUME* curVolume = (HYDRO_VOLUME*)curComp;
			
			curVolume->alpha_f = (1.0-value);
			curVolume->alpha_g = value;
			break;
		}
		case TA_BCVOIDFRACTION:
		{
			double value = lua_tonumber(L,4);
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return 0;}
			HYDRO_BOUNDARY* curVolume = (HYDRO_BOUNDARY*)curComp;
			
			curVolume->alpha_f = (1.0-value);
			curVolume->alpha_g = value;
			break;
		}
		case TA_PRESSURETEMPERATURE:
		{
			double value1 = lua_tonumber(L,4);
			double value2 = lua_tonumber(L,5);
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return 0;}
			HYDRO_VOLUME* curVolume = (HYDRO_VOLUME*)curComp;
			
			curVolume->P = value1;
			curVolume->T_f = value2;
			curVolume->T_g = value2;
			break;
		}
		case TA_BCPRESSURETEMPERATURE:
		{
			double value1 = lua_tonumber(L,4);
			double value2 = lua_tonumber(L,5);
			HYDRO_COMPONENT* curComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
			if (curComp == NULL) {return 0;}
			HYDRO_BOUNDARY* curVolume = (HYDRO_BOUNDARY*)curComp;
			
			curVolume->P = value1;
			curVolume->T_f = value2;
			curVolume->T_g = value2;
			break;
		}
		case TA_JLOSSCOEFFICIENTS:
		{
			double value1 = lua_tonumber(L,4);
			double value2 = lua_tonumber(L,5);
			HYDRO_SJUNCTION* curComp = curSystem->hydroSJunctionStack.GetItem(lcompHandle);
			
			curComp->floss = value1;
			curComp->rloss = value2;
			break;
		}
		default:
			break;
	}
	return 0;
}