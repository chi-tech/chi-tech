#include "../lbkes_k_eigenvalue_solver.h"

#include "chi_lua.h"
#include "lbkes_lua_utils.h"

#include "chi_log.h"

using namespace lbs;

//############################################################
/**Set properties for the solver.
\param SolverIndex int Handle to the solver.
\param Property string Property name for a specific property.
\param PropertyValue varying Value to set to the property

##_

###Property
"MAX_ITERATIONS" Sets the maximum k-iterations.[Default=1000]\n\n
"TOLERANCE" Sets the tolerance on k-eff for convergence.[Default=1.0e-8]\n\n

*/
int chiLBKESSetProperty(lua_State *L)
{
  const std::string fname = "chiLBKESSetProperty";
  const int num_args = lua_gettop(L);
  if (num_args < 3)
    LuaPostArgAmountError(fname, 3, num_args);

  //============================================= Get the solver
  LuaCheckNilValue(fname, L, 1);
  const int solver_handle = lua_tonumber(L, 1);

  auto& solver = chi::GetStackItem<lbs::DiscOrdKEigenvalueSolver>(chi::solver_stack,
                                                                  solver_handle,
                                                                  fname);

  //============================================= Get property name
  LuaCheckNilValue(fname, L, 2);
  const std::string property = lua_tostring(L,2);

  //============================================= Handle properties
  if (property == "MAX_ITERATIONS")
  {
    LuaCheckNilValue(fname, L, 3);
    const int max_iters = lua_tointeger(L, 3);

    if (max_iters <= 0)
    {
      chi::log.LogAllError()
          << fname << ": Invalid max_iterations value. "
          << "Must be greater than 0.";
      chi::Exit(EXIT_FAILURE);
    }
    solver.max_iterations_ = static_cast<size_t>(max_iters);

    chi::log.Log()
      << "LinearBoltzmann::KEigenvalueSolver: "
      << "max_iterations set to " << solver.max_iterations_ << ".";
  }

  else if (property == "TOLERANCE")
  {
    LuaCheckNilValue(fname, L, 3);
    const double tol = lua_tonumber(L, 3);

    if (tol < 0.0 or tol > 1.0)
    {
      chi::log.LogAllError()
          << fname << ": Invalid value for tolerance. "
          << "Must be in the range (0.0, 1.0].";
      chi::Exit(EXIT_FAILURE);
    }
    solver.tolerance_ = tol;

    char buff[100];
    snprintf(buff,100, "%.4e", tol);

    chi::log.Log()
        << "LinearBoltzmann::KEigenvalueSolver: "
        << "tolerance set to " << buff << ".";
  }
  else if (property == "K_EIGEN_METHOD")
  {
    LuaCheckNilValue(fname, L, 3);
    LuaCheckStringValue(fname, L, 3);
    const std::string method = lua_tostring(L, 3);

    solver.k_eigen_method_ = method;

    chi::log.Log()
      << "LinearBoltzmann::KEigenvalueSolver: "
      << "k_eigen_method set to " << method << ".";
  }
  else
  {
    chi::log.LogAllError() << fname << ": Invalid property index.";
    chi::Exit(EXIT_FAILURE);
  }
  return 0;
}
