#include "TimeStepController.h"

#include <cmath>

namespace chi_physics
{

chi::InputParameters TimeStepController::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddOptionalParameter("dt", 0.01, "Initial timestep to use");
  params.AddOptionalParameter("time", 0.0, "Initial time");
  params.AddOptionalParameter(
    "time_index", 0, "Time index. Useful for output control.");
  params.AddOptionalParameter("end_time", 1.0, "End time");
  params.AddOptionalParameter(
    "max_time_steps",
    -1,
    "Maximum number of timesteps to take. A negative number disables this.");

  params.AddOptionalParameter("eps",
                              1.0e-12,
                              "General time tolerance. This is used to "
                              "constrict the dt to the end-time.");

  return params;
}

TimeStepController::TimeStepController(const chi::InputParameters& params)
  : ChiObject(params),
    dt_(params.GetParamValue<double>("dt")),
    time_(params.GetParamValue<double>("time")),
    t_index_(params.GetParamValue<int>("time_index")),
    end_time_(params.GetParamValue<double>("end_time")),
    max_time_steps_(params.GetParamValue<int>("max_time_steps")),
    general_tolerance_(params.GetParamValue<double>("eps"))
{
}

/**Overridable method to get the timestep size.*/
double TimeStepController::GetTimeStepSize()
{
  dt_ = StepSizeLimitedToEndTime();

  return dt_;
}

/**Returns the current controller time.*/
double TimeStepController::Time() const
{
  return time_;
}

/**Returns the current time index.*/
int TimeStepController::TimeIndex() const
{
  return t_index_;
}

/**Provides the equivalent timestep for ending at the end-time.*/
double TimeStepController::StepSizeLimitedToEndTime() const
{
  double dt_limited = dt_;

  if ((end_time_ - time_) <= dt_) dt_limited = end_time_ - time_;

  return dt_limited;
}

/**Determines if the time is at the end time (within tolerance) or greater.*/
bool TimeStepController::AtEndTime() const
{
  if (std::fabs(time_ - end_time_) < general_tolerance_)
    return true;

  if (time_ > end_time_)
    return true;

  return false;
}

/**Manually set the time step size.*/
void TimeStepController::SetTimeStepSize(double dt) { dt_ = dt; }


/**Advances the controller's state. The most basic action here is to
  * advance time and the time index.*/
void TimeStepController::Advance()
{
  time_ += dt_;
  t_index_ += 1;

  if (AtEndTime()) finished_ = true;
}

std::string TimeStepController::StringTimeInfo() const
{
  std::stringstream outstr;
  outstr << "\nTime step " << t_index_;
  {
    char buffer[100];
    snprintf(buffer, 100, ", time = %g", time_);
    outstr << buffer;
  }
  {
    char buffer[100];
    snprintf(buffer, 100, ", dt = %g", dt_);
    outstr << buffer;
  }

  return outstr.str();
}

} // namespace chi_physics