#ifndef _chi_physicsmaterial_h
#define _chi_physicsmaterial_h

#include "../chi_physics_namespace.h"
#include <string>
#include <vector>

//Property indices
#define THERMAL_CONDUCTIVITY      0
#define SCALAR_VALUE              1
#define TRANSPORT_XSECTIONS       10
#define ISOTROPIC_MG_SOURCE       11

//###################################################################
/** Base class for material properties.*/
class chi_physics::MaterialProperty
{
public:
  int  type_index;
  std::string property_name;

  MaterialProperty()
  {
    type_index = 0;
  }
  virtual double GetScalarValue()
  {
    return 0.0;
  }
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