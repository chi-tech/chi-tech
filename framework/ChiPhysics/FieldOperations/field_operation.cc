#include "field_operation.h"

namespace chi_physics::field_operations
{

//Since there are no input parameters we will not register this object

// ##################################################################
/**Returns the input parameters.*/
chi_objects::InputParameters FieldOperation::GetInputParameters()
{
  return ChiObject::GetInputParameters();
}

// ##################################################################
/**Constructor.*/
FieldOperation::FieldOperation(const chi_objects::InputParameters& params)
  : ChiObject(params)
{
}

} // namespace chi_physics::field_operations