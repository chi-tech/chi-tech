#include "GridBasedFieldFunctionInterface.h"

#include "fieldfunction_gridbased.h"

namespace chi_physics
{

chi::InputParameters GridBasedFieldFunctionInterface::GetInputParameters()
{
  chi::InputParameters params = FieldFunctionInterface::GetInputParameters();

  return params;
}

GridBasedFieldFunctionInterface::GridBasedFieldFunctionInterface(
  const chi::InputParameters& params)
  : FieldFunctionInterface(params)
{
}

const FieldFunctionGridBased*
GridBasedFieldFunctionInterface::GetGridBasedFieldFunction() const
{
  const auto* ff_ptr = GetFieldFunction();

  const auto* grid_based_ff_ptr =
    dynamic_cast<const FieldFunctionGridBased*>(ff_ptr);

  if (not grid_based_ff_ptr)
    return nullptr;

  return grid_based_ff_ptr;
}

} // namespace chi_physics