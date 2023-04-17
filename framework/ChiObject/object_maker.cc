#include "object_maker.h"

#include "chi_runtime.h"
#include "chi_log.h"

// ###################################################################
/**Access to the singleton*/
ChiObjectMaker& ChiObjectMaker::GetInstance() noexcept
{
  static ChiObjectMaker singleton;
  return singleton;
}

// ###################################################################
/**Returns a constant reference to the object registry.*/
const std::map<std::string, ChiObjectMaker::ObjectRegistryEntry>&
ChiObjectMaker::Registry() const
{
  return object_registry_;
}

// ###################################################################
/**Makes an object with the given parameters and places on the global
 * object stack. Returns a handle to the object.*/
size_t
ChiObjectMaker::MakeObject(const chi_objects::ParameterBlock& params) const
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

  return chi::object_stack.size() - 1;
}

// ###################################################################
/**Makes an object with the given parameters and places on the global
 * object stack. Returns a handle to the object.*/
size_t
ChiObjectMaker::MakeObjectType(const std::string& type,
                               const chi_objects::ParameterBlock& params) const
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

  return chi::object_stack.size() - 1;
}

// ##################################################################
/**Dumps the registry to stdout.*/
void ChiObjectMaker::DumpRegister() const
{
  chi::log.Log() << "\n\n";
  for (const auto& [key, entry] : object_registry_)
  {
    if (chi::log.GetVerbosity() == 0)
    {
      chi::log.Log() << key;
      continue;
    }

    chi::log.Log() << "OBJECT_BEGIN " << key;

    const auto in_params = entry.get_in_params_func();
    in_params.DumpParameters();

    chi::log.Log() << "OBJECT_END\n\n";
  }
  chi::log.Log() << "\n\n";
}