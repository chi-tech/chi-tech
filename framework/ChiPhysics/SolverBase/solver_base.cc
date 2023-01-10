#include "chi_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

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