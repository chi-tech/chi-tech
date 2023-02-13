#include "lbs_common_lua_functions.h"

#define RegisterFunction(x) \
  lua_register(L, #x, lbs::common_lua_utils::x)

#define RegisterTable(x) \
        lua_newtable(L); \
        lua_setglobal(L, #x)

#define RegisterConstantValueToTable(const_name,const_value,namespace_name) \
        lua_getglobal(L,#namespace_name); \
        lua_pushstring(L,#const_name); \
        lua_pushnumber(L,const_value); \
        lua_settable(L,-3); \
        lua_pop(L,1)

namespace lbs::common_lua_utils
{
  void RegisterLuaEntities(lua_State* L)
  {
    RegisterTable(LBSGroupset);
    RegisterFunction(chiLBSCreateGroupset);
    RegisterFunction(chiLBSCreateGroup);
    RegisterFunction(chiLBSGroupsetAddGroups);
    RegisterFunction(chiLBSGroupsetSetQuadrature);
    RegisterFunction(chiLBSGroupsetSetAngleAggregationType);
    RegisterConstantValueToTable(ANGLE_AGG_SINGLE,   1,LBSGroupset);
    RegisterConstantValueToTable(ANGLE_AGG_POLAR,    2,LBSGroupset);
    RegisterConstantValueToTable(ANGLE_AGG_AZIMUTHAL,3,LBSGroupset);
    RegisterFunction(chiLBSGroupsetSetAngleAggDiv);
    RegisterFunction(chiLBSGroupsetSetGroupSubsets);
    RegisterFunction(chiLBSGroupsetSetIterativeMethod);
    RegisterFunction(chiLBSGroupsetSetResidualTolerance);
    RegisterFunction(chiLBSGroupsetSetMaxIterations);
    RegisterFunction(chiLBSGroupsetSetGMRESRestartIntvl);
    RegisterFunction(chiLBSGroupsetSetEnableSweepLog);
    RegisterFunction(chiLBSGroupsetSetWGDSA);
    RegisterFunction(chiLBSGroupsetSetTGDSA);
  }
}