#include "chi_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

chi_objects::InputParameters chi_physics::Solver::GetInputParameters()
{
  chi_objects::InputParameters params;

  params.AddRequiredParameter<std::string>("name");

  return params;
}

chi_physics::Solver::Solver(const chi_objects::InputParameters& params) :
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