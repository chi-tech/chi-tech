#include "chi_physicsmaterial.h"

//###################################################################
/** Base class method for pushing lua table.*/
void chi_physics::MaterialProperty::PushLuaTable(lua_State* L)
{
  lua_newtable(L);
  lua_pushstring(L,"is_empty");
  lua_pushboolean(L,true);
  lua_settable(L,-3);
}