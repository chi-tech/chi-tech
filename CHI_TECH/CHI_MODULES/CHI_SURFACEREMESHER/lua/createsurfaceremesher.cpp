#include "../../../CHI_LUA/chi_lua.h"
#include<math.h>

#include "../chi_surfaceremesher.h"


extern CHI_VECTOR<CHI_SURFACEREMESHER> chisurfaceRemesherStack;

//########################################################## Create empty system
/** Creates a remeshing object

\return Returns a unique handle to the remesher

\ingroup LuaSurfaceRemesher
\author Jan*/
int chiSurfaceRemesherCreate(lua_State *L)
{
	CHI_SURFACEREMESHER* newMesher = new CHI_SURFACEREMESHER;
	
	int newIndex = chisurfaceRemesherStack.PushItem(newMesher);
	
	lua_pushnumber(L,newIndex);
	return 1;
}