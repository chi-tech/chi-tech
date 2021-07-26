#include "lua_utils.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

LinearBoltzmann::KEigenvalue::Solver* LinearBoltzmann::lua_utils::
  GetKEigenvalueSolverByHandle(int handle, const std::string& calling_function_name)
{
  LinearBoltzmann::KEigenvalue::Solver* k_solver;
  try{

    k_solver = dynamic_cast<LinearBoltzmann::KEigenvalue::Solver*>(
      chi_physics_handler.solver_stack.at(handle));

    if (not k_solver)
      throw std::logic_error(calling_function_name +
      ": Invalid solver at given handle (" +
      std::to_string(handle) + "). "
      "The solver is not of type LinearBoltzmann::KEigenvalue::Solver.");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(calling_function_name + ": Invalid solver-handle (" +
                           std::to_string(handle) + ").");
  }

  return k_solver;
}