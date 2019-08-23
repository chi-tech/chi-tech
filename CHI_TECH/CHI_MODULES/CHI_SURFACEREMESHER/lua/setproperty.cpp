#include "../../../CHI_LUA/chi_lua.h"
#include<math.h>

#include "../chi_surfaceremesher.h"


extern CHI_VECTOR<CHI_SURFACEREMESHER> chisurfaceRemesherStack;

//########################################################## Create empty system
/** Creates a remeshing object

\return Returns a unique handle to the remesher

\ingroup LuaSurfaceRemesher
\author Jan*/
int chiSurfaceRemesherSetProperty(lua_State *L)
{
	int mesherIndex = lua_tonumber(L,1);
	
	CHI_SURFACEREMESHER* mesher = chisurfaceRemesherStack.GetItem(mesherIndex);
	if (mesher==NULL)
	{
		printf("Mesher not found during chiSurfaceRemesherSetProperty");
		lua_pushnumber(L,-1);
		return 1;
	}
	
	int propertyNumber = lua_tonumber(L,2);
	
	//=========================================== Property 1 MESHER_BASE_SIZE
	if (propertyNumber==1)
	{
		mesher->baseSize = lua_tonumber(L,3);
	}
	
	//=========================================== Property 2 MESHER_MIN_SIZE
	if (propertyNumber==2)
	{
		mesher->absoluteMinumumSize = lua_tonumber(L,3);
	}
	
	//=========================================== Property 3 MESHER_CONSOLIDATE_EDGES
	if (propertyNumber==3)
	{
		mesher->removeInteriorFaces = lua_toboolean(L,3);
	}
	
	//=========================================== Property 4 MESHER_KEEP2D_ORIENTATION
	if (propertyNumber==4)
	{
		mesher->keep2Dorientation = lua_toboolean(L,3);
	}

	//=========================================== Property 5 MESHER_RECALCULATE_NORMAL
	if (propertyNumber==5)
	{
		mesher->recalcNormal = lua_toboolean(L,3);
	}
	
	return 0;
}