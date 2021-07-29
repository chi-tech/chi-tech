#ifndef KEIGEN_SOLVER_H
#define KEIGEN_SOLVER_H

#include "../k_eigenvalue_solver.h"

namespace LinearBoltzmann
{
  namespace lua_utils
  {
    //###################################################################
    /** Obtains a pointer to a LinearBoltzmann::Solver object or an object
     * derived from LinearBoltzmann::Solver
     *
     * \param handle int Index in the chi_physics_handler where the solve object
     *                   should be located.
     * \param calling_function_name string The string used to print error messages,
     *                              should uniquely identify the calling function.
     *
     */
    LinearBoltzmann::KEigenvalue::Solver*
      GetKEigenvalueSolverByHandle(int handle, const std::string& calling_function_name);
  }
}

#endif