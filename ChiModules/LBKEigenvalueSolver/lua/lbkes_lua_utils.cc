#include "lbkes_lua_utils.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

LinearBoltzmann::KEigenvalueSolver* LinearBoltzmann::k_eigenvalue_lua_utils::
  GetSolverByHandle(int handle, const std::string& calling_function_name)
{
  LinearBoltzmann::KEigenvalueSolver* lbkes_solver;
  try{

    lbkes_solver = dynamic_cast<LinearBoltzmann::KEigenvalueSolver*>(
      chi_physics_handler.solver_stack.at(handle));

    if (not lbkes_solver)
      throw std::logic_error(calling_function_name +
      ": Invalid solver at given handle (" +
      std::to_string(handle) + "). "
      "The solver is not of type LinearBoltzmann::KEigenvalueSolver.");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(calling_function_name + ": Invalid solver-handle (" +
                           std::to_string(handle) + ").");
  }

  return lbkes_solver;
}
