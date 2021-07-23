#ifndef LBS_LUA_UTILS_H
#define LBS_LUA_UTILS_H

#include "../lbs_linear_boltzmann_solver.h"

namespace LinearBoltzmann
{
  namespace lua_utils
  {
    //###################################################################
    /** Obtains a pointer to a LinearBoltzmann::Solver object or an object
     * derived from LinearBoltzmann::Solver
     *
     * \param handle int index in the chi_physics_handler where the solve object should be located
     * \param calling_function_name string used to when print error messages, should uniquely
     *                              identify the calling function
     *
     */
    LinearBoltzmann::Solver* GetSolverByHandle(int handle, const std::string& calling_function_name);
  }
}

#endif
