#include "../../../CHI_LUA/chi_lua.h"
#include<math.h>

#include "../chi_surfaceremesher.h"

extern CHI_VECTOR<CHI_SURFACEREMESHER> chisurfaceRemesherStack;
extern CHI_VECTOR<CHI_SURFACE> chisurfaceMeshStack;

/** \defgroup LuaSurfaceRemesher Surface Remesher
 *\ingroup LuaGeneralUtilities*/

//########################################################## Create empty system
/** Remeshes the given surface

\return Returns a unique handle to a new surface

\ingroup LuaSurfaceRemesher
\author Jan*/
int chiSurfaceRemesherExecuteMeshing(lua_State *L)
{
    int numArgs = lua_gettop(L);
    
    int mesherIndex = lua_tonumber(L,1);
    
    CHI_SURFACEREMESHER* mesher = chisurfaceRemesherStack.GetItem(mesherIndex);
    if (mesher==NULL)
    {
        printf("Mesher not found during chiSurfaceRemesherExecuteMeshing");
        lua_pushnumber(L,-1);
        return 1;
    }
    
    int surfaceIndex = lua_tonumber(L,2);
    
    CHI_SURFACE* curSurface = chisurfaceMeshStack.GetItem(surfaceIndex);
    
    if (curSurface==NULL)
    {
        printf("Surface not found during chiSurfaceRemesherExecuteMeshing");
        lua_pushnumber(L,-1);
        return 1;
    }
    
    
    if (numArgs>=3)
    {
        mesher->removeInteriorFaces = lua_toboolean(L,2);
    }
    if (numArgs==4)
    {
        mesher->keep2Dorientation = lua_toboolean(L,3);
    }

    mesher->ExecuteMeshing(curSurface);
    
    int newSurf = chisurfaceMeshStack.PushItem(mesher->remeshedSurface);
    
    lua_pushnumber(L,newSurf);

    return 1;
}
