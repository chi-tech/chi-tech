#include "../../../ChiLua/chi_lua.h"
#include "../chi_thermoalpha.h"
#include<math.h>

#include "../hydro_volume/hydro_volume.h"
#include "../hydro_boundary/hydro_boundary.h"

extern CHI_VECTOR<CHI_THERMOSYSTEM> chithermoSystemStack;

/** \defgroup LuaThermoalpha Thermoalpha
 * \ingroup LuaModules*/

//########################################################## Create empty system
/** void chiThermoCreateSystem() Generate an empty Thermal-hydraulic system.

\return Returns a unique handle for the created system.

\ingroup LuaThermoalpha
\author Jan*/
int chiThermoCreateSystem(lua_State *L)
{
	CHI_THERMOSYSTEM* newSystem = new CHI_THERMOSYSTEM;
	
	int index = chithermoSystemStack.PushItem(newSystem);
	
	lua_pushnumber(L,index);
	
	return 1;
}










//######################################################### Connect Components
/** void chiThermoConnectTwoComponents() Connect two hydrodynamic components using a single junction.

\param systemHandle  Handle to the system to which all the components belong.
\param leftComponent Component to be connected to the left of the junction.
\param sjunc Single junction to be used for the connection.
\param rigtComponent Component to be connected to the right of the junction.
\param mode 0=end-begin, 1=begin-end, 2=end-end,

\ingroup LuaThermoalpha
\author Jan*/
int chiThermoConnectTwoComponents(lua_State *L)
{
	//===================================================== Getting the system
	int systemHandle = lua_tonumber(L,1);
	CHI_THERMOSYSTEM* curSystem = chithermoSystemStack.GetItem(systemHandle);
	
	if (curSystem == NULL) {return 0;}
	
	//===================================================== Getting the left component
	int lcompHandle = lua_tonumber(L,2);
	HYDRO_COMPONENT* lComp = curSystem->hydroComponentStack.GetItem(lcompHandle);
	
	if (lComp == NULL) {return 0;}
	
	//===================================================== Getting the junction
	int juncHandle = lua_tonumber(L,3);
	HYDRO_SJUNCTION* junc = curSystem->hydroSJunctionStack.GetItem(juncHandle);
	
	if (junc == NULL) {return 0;}
	
	//===================================================== Getting the right component
	int rcompHandle = lua_tonumber(L,4);
	HYDRO_COMPONENT* rComp = curSystem->hydroComponentStack.GetItem(rcompHandle);
	
	if (rComp == NULL) {return 0;}
	
	int mode = lua_tonumber(L, 5);

	if (lComp->componentType==1)
	{
		HYDRO_VOLUME* realHydro = (HYDRO_VOLUME*)lComp;
		realHydro->rigtJunctionIndex = juncHandle;
	}
	else
	{
		HYDRO_BOUNDARY* realHydro = (HYDRO_BOUNDARY*)lComp;
		realHydro->junctionIndex = juncHandle;
	}
	
	if (rComp->componentType==1)
	{
		HYDRO_VOLUME* realHydro = (HYDRO_VOLUME*)rComp;
		realHydro->leftJunctionIndex = juncHandle;
	}
	else
	{
		HYDRO_BOUNDARY* realHydro = (HYDRO_BOUNDARY*)rComp;
		realHydro->junctionIndex = juncHandle;
	}
	
	junc->leftHComponent = lComp;
	junc->rigtHComponent = rComp;
	
	
	return 0;
}



//######################################################### Initialize System
/** bool chiThermoInitialize() Initializes system.

\param systemHandle int Handle to the system which should be initialized.

\return Returns true if successfully initialized and false otherwise.

\ingroup LuaThermoalpha
\author Jan*/
int chiThermoInitialize(lua_State *L)
{
	//===================================================== Getting the system
	int systemHandle = lua_tonumber(L,1);
	CHI_THERMOSYSTEM* curSystem = chithermoSystemStack.GetItem(systemHandle);
	
	
	lua_pushboolean(L,curSystem->Initialize());
	return 1;
}
