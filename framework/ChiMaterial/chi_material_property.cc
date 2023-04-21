#include "chi_material_property.h"

#include "ChiObject/object_maker.h"

namespace chi_objects
{

RegisterChiObject(chi_objects, MaterialProperty);

InputParameters MaterialProperty::GetInputParameters()
{
  InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>(
    "name", "Text name associated with this property");

  return params;
}

MaterialProperty::MaterialProperty(const chi_objects::InputParameters& params)
  : name_(params.GetParamValue<std::string>("name"))
{
}

const std::string& MaterialProperty::TextName() const { return name_; }

} // namespace chi_objects