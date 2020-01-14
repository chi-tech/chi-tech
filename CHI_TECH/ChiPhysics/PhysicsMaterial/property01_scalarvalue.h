#ifndef _chi_physics_property_scalarvalue_h
#define _chi_physics_property_scalarvalue_h

#include "chi_physicsmaterial.h"

//###################################################################
/**Simple scalar material property.*/
class chi_physics::ScalarValue : public chi_physics::MaterialProperty
{
public:
  double value=1.0;

  ScalarValue() : MaterialProperty(PropertyType::SCALAR_VALUE) {}

  double GetScalarValue() {return value;}
  void PushLuaTable(lua_State* L) override
  {
    lua_newtable(L);
    lua_pushstring(L,"is_empty");
    lua_pushboolean(L,false);
    lua_settable(L,-3);

    lua_pushstring(L,"value");
    lua_pushnumber(L,value);
    lua_settable(L,-3);
  }

};

#endif