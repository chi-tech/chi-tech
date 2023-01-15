#include "ffinterpol_lua.h"

#include "chi_lua.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)
#define LUA_CMACRO1(x,y) \
        lua_pushnumber(L, y); \
        lua_setglobal(L, #x)

void chi_mesh::ff_interpolation_lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiFFInterpolationCreate);
    LUA_CMACRO1(SLICE , 1);
    LUA_CMACRO1(LINE  , 2);
    LUA_CMACRO1(VOLUME, 3);
  LUA_FMACRO1(chiFFInterpolationSetProperty);
    LUA_CMACRO1(SLICE_POINT,   1);
    LUA_CMACRO1(SLICE_NORMAL,   2);
    LUA_CMACRO1(SLICE_TANGENT,   3);
    LUA_CMACRO1(SLICE_BINORM,   4);
    LUA_CMACRO1(OPERATION,   5);
      LUA_CMACRO1(OP_SUM,   10);
      LUA_CMACRO1(OP_AVG,   11);
      LUA_CMACRO1(OP_MAX,   12);
      LUA_CMACRO1(OP_SUM_LUA,   13);
      LUA_CMACRO1(OP_AVG_LUA,   14);
      LUA_CMACRO1(OP_MAX_LUA,   15);
    LUA_CMACRO1(LOGICAL_VOLUME,   8);

    LUA_CMACRO1(ADD_FIELDFUNCTION,   9);
    LUA_CMACRO1(SET_FIELDFUNCTIONS,   10);

    LUA_CMACRO1(LINE_FIRSTPOINT,   11);
    LUA_CMACRO1(LINE_SECONDPOINT,   12);
    LUA_CMACRO1(LINE_NUMBEROFPOINTS,   13);
    LUA_CMACRO1(LINE_CUSTOM_ARRAY,   14);
  LUA_FMACRO1(chiFFInterpolationInitialize);
  LUA_FMACRO1(chiFFInterpolationExecute);
  LUA_FMACRO1(chiFFInterpolationExportPython);
  LUA_FMACRO1(chiFFInterpolationGetValue);
}