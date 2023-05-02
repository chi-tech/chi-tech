#include "chi_object.h"

/**Returns the input parameters.*/
chi_objects::InputParameters ChiObject::GetInputParameters()
{
  return {}; // Returns an empty block
}

ChiObject::ChiObject() {}

ChiObject::ChiObject(const chi_objects::InputParameters&)
{
}

void ChiObject::SetStackID(size_t stack_id) { stack_id_ = stack_id; }

void ChiObject::SetParamBlockUsedAtConstruction(
  const chi_objects::ParameterBlock& params)
{
  param_block_used_at_construction_ = params;
}

size_t ChiObject::StackID() const { return stack_id_; }

const chi_objects::ParameterBlock&
ChiObject::ParamBlockUsedAtConstruction() const
{
  return param_block_used_at_construction_;
}

void ChiObject::PushOntoStack(std::shared_ptr<ChiObject>& new_object)
{
  chi::object_stack.push_back(new_object);
  new_object->SetStackID(chi::object_stack.size() - 1);
}