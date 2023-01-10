#include "lbs_curvilinear_solver_lua.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)

#define LUA_CTABLE1(x) \
        lua_newtable(L); \
        lua_setglobal(L, #x)

#define LUA_CADDCONST_VALUE_TO_TABLE1(const_name,const_value,namespace_name) \
        lua_getglobal(L,#namespace_name); \
        lua_pushstring(L,#const_name); \
        lua_pushnumber(L,const_value); \
        lua_settable(L,-3); \
        lua_pop(L,1)

void lbs_curvilinear::lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiLBSCurvilinearCreateSolver);

  LUA_CTABLE1(LBSCurvilinear);
  LUA_CADDCONST_VALUE_TO_TABLE1(CYLINDRICAL, 2, LBSCurvilinear);
  LUA_CADDCONST_VALUE_TO_TABLE1(SPHERICAL,   3, LBSCurvilinear);
}