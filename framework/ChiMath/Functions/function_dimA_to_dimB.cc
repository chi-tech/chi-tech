#include "function_dimA_to_dimB.h"

namespace chi_math
{

chi_objects::InputParameters FunctionDimAToDimB::GetInputParameters()
{
  chi_objects::InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<size_t>(
    "input_dimension",
    "The dimension of the input values (excluding the position).");

  params.AddRequiredParameter<size_t>(
    "output_dimension",
    "The dimension of the output values (excluding the position).");

  return params;
}

FunctionDimAToDimB::FunctionDimAToDimB(
  const chi_objects::InputParameters& params) :
  ChiObject(params),
  input_dimension_(params.GetParamValue<size_t>("input_dimension")),
  output_dimension_(params.GetParamValue<size_t>("output_dimension"))
{
}

} // namespace chi_math