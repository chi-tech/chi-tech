#ifndef _chi_physics_property_thermalconductivity_h
#define _chi_physics_property_thermalconductivity_h

#include "chi_physicsmaterial.h"

//###################################################################
/** Basic thermal conductivity material property.*/
class chi_physics::ThermalConductivity : public chi_physics::MaterialProperty
{
public:
  double k;

  ThermalConductivity()
  {
    type_index = THERMAL_CONDUCTIVITY;
    k = 1.0;
  }
  double GetScalarValue()
  {
    return k;
  }
};



#endif