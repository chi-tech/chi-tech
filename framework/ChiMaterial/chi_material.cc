#include "chi_material.h"

#include "ChiObject/object_maker.h"

namespace chi_objects
{

RegisterChiObject(chi_objects, Material);

InputParameters Material::GetInputParameters()
{
  InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>(
    "name", "The text name that will be associated with this material.");

  params.AddRequiredParameterArray(
    "properties",
    "Expects an array object handles that represents the properties.");

  return params;
}

Material::Material(const chi_objects::InputParameters& params)
  : name_(params.GetParamValue<std::string>("name"))
{
}

} // namespace chi_objects