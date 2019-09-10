#ifndef _chi_physics_property_scalarvalue_h
#define _chi_physics_property_scalarvalue_h

#include "chi_physicsmaterial.h"

//###################################################################
/**Simple scalar material property.*/
class chi_physics::ScalarValue : public chi_physics::MaterialProperty
{
public:
  double value;

  ScalarValue()
  {
    type_index = SCALAR_VALUE;
    value = 1.0;
  }

  double GetScalarValue()
  {
    return value;
  }
};

#endif