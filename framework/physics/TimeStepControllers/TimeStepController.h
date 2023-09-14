#ifndef CHITECH_TIMESTEPCONTROLLER_H
#define CHITECH_TIMESTEPCONTROLLER_H

#include "ChiObject.h"

namespace chi_physics
{

/**Base class for all timestep controllers.*/
class TimeStepController : public ChiObject
{
public:
  virtual double GetTimeStepSize();
  void SetTimeStepSize(double dt);

protected:
  static chi::InputParameters GetInputParameters();
  explicit TimeStepController(const chi::InputParameters& params);

  double current_timestep_size_;
};

}

#endif // CHITECH_TIMESTEPCONTROLLER_H
