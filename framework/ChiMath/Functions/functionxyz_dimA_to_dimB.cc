#include "functionxyz_dimA_to_dimB.h"

#include "ChiObject/object_maker.h"

namespace chi_math::functions
{

RegisterChiObject(chi_math::functions, FunctionXYZDimAToDimB);

chi_objects::InputParameters FunctionXYZDimAToDimB::GetInputParameters()
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

FunctionXYZDimAToDimB::FunctionXYZDimAToDimB(
  const chi_objects::InputParameters& params)
  : ChiObject(params),
    input_dimension_(params.GetParamValue<size_t>("input_dimension")),
    output_dimension_(params.GetParamValue<size_t>("output_dimension"))
{
}

size_t FunctionXYZDimAToDimB::InputDimension() const { return input_dimension_; }

size_t FunctionXYZDimAToDimB::OutputDimension() const
{
  return output_dimension_;
}

std::vector<double> FunctionXYZDimAToDimB::Evaluate(
  double x, double y, double z, const std::vector<double>& values) const
{
  ChiLogicalError(" method not implemented");
  return {};
}

} // namespace chi_math::functions