#ifndef CHI_PHYSICS_PROPERTY_SCALAR_VALUE_H
#define CHI_PHYSICS_PROPERTY_SCALAR_VALUE_H

#include "material_property_base.h"


namespace chi_physics
{

//###################################################################
/**Simple scalar material property.*/
class ScalarValue : public chi_physics::MaterialProperty
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

}//namespace chi_physics

#endif //CHI_PHYSICS_PROPERTY_SCALAR_VALUE_H