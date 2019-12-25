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

};

#endif