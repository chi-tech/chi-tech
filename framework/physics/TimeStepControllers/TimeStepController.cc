#include "TimeStepController.h"

namespace chi_physics
{

chi::InputParameters TimeStepController::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddOptionalParameter("initial_dt", 0.01, "Initial timestep to use");

  return params;
}

TimeStepController::TimeStepController(const chi::InputParameters& params)
  : ChiObject(params),
    current_timestep_size_(params.GetParamValue<double>("initial_dt"))
{
}

} // namespace chi_physics