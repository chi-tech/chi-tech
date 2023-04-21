#ifndef CHITECH_CHI_MATERIAL_PROPERTY_H
#define CHITECH_CHI_MATERIAL_PROPERTY_H

#include "ChiObject/chi_object.h"

namespace chi_objects
{

/**Base class for a material property.*/
class MaterialProperty : public ChiObject
{
private:
  const std::string name_;

public:
  static chi_objects::InputParameters GetInputParameters();
  explicit MaterialProperty(const chi_objects::InputParameters& params);

  const std::string& TextName() const;
};

}

#endif // CHITECH_CHI_MATERIAL_PROPERTY_H
