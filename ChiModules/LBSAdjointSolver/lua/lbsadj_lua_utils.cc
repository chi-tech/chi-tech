#include "lbsadj_lua_utils.h"

#include "chi_runtime.h"

//###################################################################
/** Obtains a pointer to a lbs_adjoint::AdjointSolver object or an object
 * derived from lbs_adjoint::AdjointSolver
 *
 * \param handle int Index in the chi_physics_handler where the solve object
 *                   should be located.
 * \param calling_function_name string The string used to print error messages,
 *                              should uniquely identify the calling function.
 *
 */
lbs_adjoint::AdjointSolver& lbs_adjoint::lua_utils::GetSolverByHandle(
  int handle, const std::string& calling_function_name)
{
  std::shared_ptr<lbs_adjoint::AdjointSolver> adj_solver;
  try{

    adj_solver = std::dynamic_pointer_cast<lbs_adjoint::AdjointSolver>(
      chi::solver_stack.at(handle));

    if (not adj_solver)
      throw std::logic_error(calling_function_name +
      ": Invalid solver at given handle (" +
      std::to_string(handle) + "). "
                               "The solver is not of type LinearBoltzmann::Solver.");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(calling_function_name + ": Invalid solver-handle (" +
    std::to_string(handle) + ").");
  }

  return *adj_solver;
}

void lbs_adjoint::lua_utils::RegisterLuaEntities(lua_State* L)
{
  lua_register(L, "chiAdjointSolverCreate",
               lbs_adjoint::lua_utils::chiAdjointSolverCreate);
  lua_register(L, "chiAdjointSolverAddResponseFunction",
               lbs_adjoint::lua_utils::chiAdjointSolverAddResponseFunction);
  lua_register(L, "chiAdjointSolverMakeExpRepFromP1Moments",
               lbs_adjoint::lua_utils::chiAdjointSolverMakeExpRepFromP1Moments);
  lua_register(L, "chiAdjointSolverExportImportanceMapBinary",
               lbs_adjoint::lua_utils::chiAdjointSolverExportImportanceMapBinary);
  lua_register(L, "chiAdjointSolverComputeInnerProduct",
               lbs_adjoint::lua_utils::chiAdjointSolverComputeInnerProduct);
  lua_register(L, "chiAdjointSolverReadFluxMomentsToBuffer",
               lbs_adjoint::lua_utils::chiAdjointSolverReadFluxMomentsToBuffer);
  lua_register(L, "chiAdjointSolverApplyFluxMomentBuffer",
               lbs_adjoint::lua_utils::chiAdjointSolverApplyFluxMomentBuffer);
}