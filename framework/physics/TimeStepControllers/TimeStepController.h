#ifndef CHITECH_TIMESTEPCONTROLLER_H
#define CHITECH_TIMESTEPCONTROLLER_H

#include "ChiObject.h"

namespace chi_physics
{

enum class TimeStepStatus
{
  SUCCESS = 0,
  FAILURE = 1,
  NEUTRAL = 2
};

/**Base class for all timestep controllers.*/
class TimeStepController : public ChiObject
{
public:
  /**Overridable method to get the timestep size.*/
  virtual double GetTimeStepSize();
  /**Returns the current controller time.*/
  double Time() const;
  /**Returns the current time index.*/
  int TimeIndex() const;

  /**Manually set the time step size.*/
  void SetTimeStepSize(double dt);

  /**Advances the controller's state. The most basic action here is to
  * advance time and the time index.*/
  virtual void Advance();

  /**Adapts according to the timestep status.*/
  virtual bool Adapt(TimeStepStatus time_step_status) {return false;}

  /**Builds a formatted string of the time information.*/
  std::string StringTimeInfo() const;

protected:
  static chi::InputParameters GetInputParameters();
  explicit TimeStepController(const chi::InputParameters& params);

  /**Provides the equivalent timestep for ending at the end-time.*/
  double StepSizeLimitedToEndTime() const;

  /**Determines if the time is at the end time (within tolerance) or greater.*/
  bool AtEndTime() const;

  // Input parameters
  double dt_;
  double time_;
  int t_index_;

  double end_time_;
  int max_time_steps_;

  const double general_tolerance_;

  // Runtime  variables
  bool finished_ = false;
};

}

#endif // CHITECH_TIMESTEPCONTROLLER_H
