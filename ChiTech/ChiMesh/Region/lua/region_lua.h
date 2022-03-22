#ifndef CHITECH_REGION_LUA_H
#define CHITECH_REGION_LUA_H

#include "chi_lua.h"

int chiRegionCreate(lua_State* L);
int chiRegionExportMeshToPython(lua_State* L);
int chiRegionExportMeshToObj(lua_State* L);
int chiRegionExportMeshToVTK(lua_State* L);

#endif //CHITECH_REGION_LUA_H
