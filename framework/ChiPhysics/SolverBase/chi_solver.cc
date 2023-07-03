#include "chi_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

/**Returns the input parameters.*/
chi::InputParameters chi_physics::Solver::GetInputParameters()
{
  chi::InputParameters params =
    ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>(
    "name",
    "A text name to associate with the solver. This name will be used "
    "in status messages and verbose iterative convergence monitors.");

  return params;
}

chi_physics::Solver::Solver(const chi::InputParameters& params)
  : ChiObject(params),
    text_name_(params.GetParamValue<std::string>("name"))
{
}

void chi_physics::Solver::Initialize()
{
  Chi::log.Log() << "\"Initialize()\" method not defined for " << TextName();
}

void chi_physics::Solver::Execute()
{
  Chi::log.Log() << "\"Execute()\" method not defined for " << TextName();
}

void chi_physics::Solver::Step()
{
  Chi::log.Log() << "\"Step()\" method not defined for " << TextName();
}

void chi_physics::Solver::Advance()
{
  Chi::log.Log() << "\"Advance()\" method not defined for " << TextName();
}