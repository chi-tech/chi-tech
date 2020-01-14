#ifndef _chi_physics_property_isotropic_mg_src_h
#define _chi_physics_property_isotropic_mg_src_h

#include "chi_physicsmaterial.h"

//###################################################################
/** Basic thermal conductivity material property.*/
class chi_physics::IsotropicMultiGrpSource : public chi_physics::MaterialProperty
{
public:
  std::vector<double> source_value_g;

  IsotropicMultiGrpSource() :
    MaterialProperty(PropertyType::ISOTROPIC_MG_SOURCE) {}

  void PushLuaTable(lua_State* L) override
  {
    lua_newtable(L);
    lua_pushstring(L,"is_empty");
    lua_pushboolean(L,false);
    lua_settable(L,-3);

    lua_pushstring(L,"G");
    lua_pushnumber(L,source_value_g.size());
    lua_settable(L,-3);

    lua_pushstring(L,"source_value_g");
    lua_newtable(L);
    int g=0;
    for (auto val : source_value_g)
    {
      ++g;
      lua_pushnumber(L,g);
      lua_pushnumber(L,val);
      lua_settable(L,-3);
    }
    lua_settable(L,-3);
  }
};

#endif