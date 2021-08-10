#include "physics_lua_utils.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

chi_physics::Solver* chi_physics::lua_utils::
  GetSolverByHandle(int handle, const std::string &calling_function_name)
{
  chi_physics::Solver* solver;
  try{

    solver = dynamic_cast<chi_physics::Solver*>(
      chi_physics_handler.solver_stack.at(handle));

    if (not solver)
      throw std::logic_error(calling_function_name +
      ": Invalid solver at given handle (" +
      std::to_string(handle) + "). "
      "The solver is not of type chi_physics::Solver.");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(calling_function_name + ": Invalid solver-handle (" +
                           std::to_string(handle) + ").");
  }

  return solver;
}