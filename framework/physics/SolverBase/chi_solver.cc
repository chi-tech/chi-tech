#include "chi_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "physics/TimeStepControllers/ConstantTimeStepController.h"

#include "ChiObjectFactory.h"

namespace chi_physics
{

/**Returns the input parameters.*/
chi::InputParameters Solver::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>(
    "name",
    "A text name to associate with the solver. This name will be used "
    "in status messages and verbose iterative convergence monitors.");

  params.AddOptionalParameter("dt", 0.01, "Desired initial timestep size.");
  params.AddOptionalParameter("time", 0.0, "Current time of the solver.");
  params.AddOptionalParameter(
    "end_time", 1.0, "Transient end-time if applicable.");
  params.AddOptionalParameter(
    "max_time_steps",
    -1,
    "Maximum number of timesteps to allow. Negative values disables this.");

  params.AddOptionalParameter(
    "timestep_controller", 0, "Timestep controller to use for timestepping.");

  using namespace chi_data_types;
  params.ConstrainParameterRange(
    "dt", AllowableRangeLowHighLimit::New(1.0e-12, 100.0));

  return params;
}

Solver::Solver(const chi::InputParameters& params)
  : ChiObject(params), text_name_(params.GetParamValue<std::string>("name")),
    dt_(params.GetParamValue<double>("dt")),
    time_(params.GetParamValue<double>("time")),
    end_time_(params.GetParamValue<double>("end_time")),
    max_time_steps_(params.GetParamValue<int>("max_time_steps"))
{
  const auto& user_params = params.ParametersAtAssignment();

  if (not user_params.Has("timestep_controller"))
  {
    chi::ParameterBlock ts_params;
    ts_params.AddParameter("initial_dt", dt_);

    size_t handle = ChiObjectFactory::GetInstance().MakeRegisteredObjectOfType(
      "chi_physics::ConstantTimeStepController", ts_params);
    time_step_controller_ = Chi::GetStackItemPtrAsType<TimeStepController>(
      Chi::object_stack, handle, __FUNCTION__);
  }
  else
  {
    const size_t handle = params.GetParamValue<size_t>("timestep_controller");
    time_step_controller_ = Chi::GetStackItemPtrAsType<TimeStepController>(
      Chi::object_stack, handle, __FUNCTION__);
  }
}

std::string Solver::TextName() const { return text_name_; }

BasicOptions& Solver::GetBasicOptions() { return basic_options_; }

const BasicOptions& Solver::GetBasicOptions() const { return basic_options_; }

std::vector<std::shared_ptr<FieldFunctionGridBased>>&
Solver::GetFieldFunctions()
{
  return field_functions_;
}

double Solver::DeltaT() const { return dt_; }
double Solver::Time() const { return time_; }
double Solver::EndTime() const { return end_time_; }

int Solver::MaxTimeSteps() const { return max_time_steps_; }

size_t Solver::TimeStepIndex() const { return t_index_; }

const std::vector<std::shared_ptr<FieldFunctionGridBased>>&
Solver::GetFieldFunctions() const
{
  return field_functions_;
}

void Solver::Initialize()
{
  Chi::log.Log() << "\"Initialize()\" method not defined for " << TextName();
}

void Solver::Execute()
{
  Chi::log.Log() << "\"Execute()\" method not defined for " << TextName();
}

void Solver::Step()
{
  Chi::log.Log() << "\"Step()\" method not defined for " << TextName();
}

void Solver::Advance()
{
  Chi::log.Log() << "\"Advance()\" method not defined for " << TextName();
}

chi::ParameterBlock Solver::GetInfo(const chi::ParameterBlock& params) const
{
  return chi::ParameterBlock{};
}

chi::ParameterBlock
Solver::GetInfoWithPreCheck(const chi::ParameterBlock& params) const
{
  if (not params.Has("name"))
  {
    Chi::log.LogAllWarning() << "chi_physics::Solver::GetInfo called without "
                                "\"name\" in the parameter list";
    return chi::ParameterBlock{};
  }
  return GetInfo(params);
}

} // namespace chi_physics
