#include "object_maker.h"

#include "chi_runtime.h"

namespace chi_objects
{

// ###################################################################
/**Access to the singleton*/
ObjectMaker& ObjectMaker::GetInstance() noexcept
{
  static ObjectMaker singleton;
  return singleton;
}

// ###################################################################
/**Returns a constant reference to the object registry.*/
const std::map<std::string, ObjectMaker::ObjectRegistryEntry>&
ObjectMaker::Registry() const
{
  return object_registry_;
}

// ###################################################################
/**Makes an object with the given parameters and places on the global
 * object stack. Returns a handle to the object.*/
size_t ObjectMaker::MakeObject(const ParameterBlock& params) const
{
  const std::string fname = __PRETTY_FUNCTION__;
  if (not params.Has("type"))
    throw std::invalid_argument(fname + ": Requires a parameter block with "
                                        "a field called \"type\".");

  const auto type = params.GetParamValue<std::string>("type");

  if (object_registry_.count(type) == 0)
    throw std::logic_error(fname + ": No registered type \"" + type +
                           "\" found.");

  auto object_entry = object_registry_.at(type);

  auto input_params = object_entry.get_in_params_func();

  input_params.SetObjectType(type);
  input_params.AssignParameters(params);

  auto new_object = object_entry.constructor_func(input_params);

  chi::object_stack.push_back(new_object);

  return chi::object_stack.size()-1;
}

// ###################################################################
/**Makes an object with the given parameters and places on the global
 * object stack. Returns a handle to the object.*/
size_t ObjectMaker::MakeObjectType(const std::string& type,
                                   const ParameterBlock& params) const
{
  const std::string fname = __PRETTY_FUNCTION__;

  if (object_registry_.count(type) == 0)
    throw std::logic_error(fname + ": No registered type \"" + type +
                           "\" found.");

  auto object_entry = object_registry_.at(type);

  auto input_params = object_entry.get_in_params_func();

  input_params.SetObjectType(type);
  input_params.AssignParameters(params);

  auto new_object = object_entry.constructor_func(input_params);

  chi::object_stack.push_back(new_object);

  return chi::object_stack.size()-1;
}

} // namespace chi_objects