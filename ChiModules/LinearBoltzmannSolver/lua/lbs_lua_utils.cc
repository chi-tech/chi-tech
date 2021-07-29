#include "lbs_lua_utils.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

LinearBoltzmann::Solver* LinearBoltzmann::lua_utils::
  GetSolverByHandle(int handle, const std::string& calling_function_name)
{
  LinearBoltzmann::Solver* lbs_solver;
  try{

    lbs_solver = dynamic_cast<LinearBoltzmann::Solver*>(
      chi_physics_handler.solver_stack.at(handle));

    if (not lbs_solver)
      throw std::logic_error(calling_function_name +
      ": Invalid solver at given handle (" +
      std::to_string(handle) + "). "
      "The solver is not of type LinearBoltzmann::Solver.");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(calling_function_name + ": Invalid solver-handle (" +
                           std::to_string(handle) + ").");
  }

  return lbs_solver;
}
