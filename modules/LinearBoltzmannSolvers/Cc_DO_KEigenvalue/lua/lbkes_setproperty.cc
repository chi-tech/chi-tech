#include "../lbkes_k_eigenvalue_solver.h"

#include "chi_lua.h"
#include "lbkes_lua_utils.h"

#include "chi_runtime.h"
#include "chi_log.h"

using namespace lbs;

//############################################################
/**Deprecated use basic options.
*/
int chiLBKESSetProperty(lua_State *L)
{
  const std::string fname = "chiLBKESSetProperty";
  chi::log.LogAllError() << fname << ": Has been deprecated. Use"
                                     "chiSolverSetBasicOption.";
  chi::Exit(EXIT_FAILURE);

  return 0;
}
