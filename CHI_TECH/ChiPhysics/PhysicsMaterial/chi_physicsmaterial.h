#ifndef _chi_physicsmaterial_h
#define _chi_physicsmaterial_h

#include "../chi_physics_namespace.h"
#include <string>
#include <vector>

#include "ChiLua/chi_lua.h"

namespace chi_physics
{
  enum class PropertyType
  {
    SCALAR_VALUE = 1,
    TRANSPORT_XSECTIONS = 10,
    ISOTROPIC_MG_SOURCE = 11
  };
}

//###################################################################
/** Base class for material properties.*/
class chi_physics::MaterialProperty
{
private:
  const PropertyType type;
public:
  std::string property_name;

  explicit MaterialProperty(PropertyType in_type) : type(in_type) {}

  virtual ~MaterialProperty() {}

  PropertyType Type() {return type;}

  virtual double GetScalarValue() { return 0.0; }
  virtual void PushLuaTable(lua_State* L);
};

//###################################################################
/** Base class for materials used in physics simulations.*/
class chi_physics::Material
{
public:
  std::vector<MaterialProperty*> properties;
public:
  std::string name;

};










#endif