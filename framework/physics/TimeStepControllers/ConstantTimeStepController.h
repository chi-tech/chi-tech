#ifndef CHITECH_CONSTANTTIMESTEPCONTROLLER_H
#define CHITECH_CONSTANTTIMESTEPCONTROLLER_H

#include "TimeStepController.h"

namespace chi_physics
{

/**Timestep controller that does not dynamically change.*/
class ConstantTimeStepController : public TimeStepController
{
public:
  static chi::InputParameters GetInputParameters();
  explicit ConstantTimeStepController(const chi::InputParameters& params);
};

}

#endif // CHITECH_CONSTANTTIMESTEPCONTROLLER_H
