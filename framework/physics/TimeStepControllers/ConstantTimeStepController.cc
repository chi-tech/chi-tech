#include "ConstantTimeStepController.h"

#include "ChiObjectFactory.h"

namespace chi_physics
{

RegisterChiObject(chi_physics, ConstantTimeStepController);

chi::InputParameters ConstantTimeStepController::GetInputParameters()
{
  chi::InputParameters params = TimeStepController::GetInputParameters();

  params.SetGeneralDescription(
    "Timestep controller that does not dynamically change.");
  params.SetDocGroup("doc_TimeStepControllers");

  return params;
}

ConstantTimeStepController::ConstantTimeStepController(
  const chi::InputParameters& params)
  : TimeStepController(params)
{
}

} // namespace chi_physics