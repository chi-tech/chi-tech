#include "lbs_disc_ord_ss_lua_utils.h"

#define RegisterFunction(x) \
  lua_register(L, #x, x)

void lbs::disc_ord_steady_state_lua_utils::RegisterLuaEntities(lua_State *L)
{
  RegisterFunction(chiLBSCreateSolver);

  RegisterFunction(chiLBSComputeBalance);
}
