#include "function_dimA_to_dimB.h"

namespace chi_math
{

chi::InputParameters FunctionDimAToDimB::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<size_t>(
    "input_dimension",
    "The dimension of the input values (excluding the position).");

  params.AddRequiredParameter<size_t>(
    "output_dimension",
    "The dimension of the output values (excluding the position).");

  return params;
}

FunctionDimAToDimB::FunctionDimAToDimB(
  const chi::InputParameters& params) :
  ChiObject(params),
  input_dimension_(params.GetParamValue<size_t>("input_dimension")),
  output_dimension_(params.GetParamValue<size_t>("output_dimension"))
{
}

} // namespace chi_math