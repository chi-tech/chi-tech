#include "chi_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

/**Returns the input parameters.*/
chi_objects::InputParameters chi_physics::Solver::GetInputParameters()
{
  chi_objects::InputParameters params =
    ChiObject::GetInputParameters();

  params.AddRequiredParameter<std::string>(
    "name",
    "A text name to associate with the solver. This name will be used "
    "in status messages and verbose iterative convergence monitors.");

  return params;
}

chi_physics::Solver::Solver(const chi_objects::InputParameters& params)
  : ChiObject(params),
    text_name_(params.GetParamValue<std::string>("name"))
{
}

void chi_physics::Solver::Initialize()
{
  chi::log.Log() << "\"Initialize()\" method not defined for " << TextName();
}

void chi_physics::Solver::Execute()
{
  chi::log.Log() << "\"Execute()\" method not defined for " << TextName();
}

void chi_physics::Solver::Step()
{
  chi::log.Log() << "\"Step()\" method not defined for " << TextName();
}

void chi_physics::Solver::Advance()
{
  chi::log.Log() << "\"Advance()\" method not defined for " << TextName();
}