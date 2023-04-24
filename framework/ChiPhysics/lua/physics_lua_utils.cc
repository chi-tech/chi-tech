#include "physics_lua_utils.h"

#include "ChiPhysics/FieldFunction/lua/fieldfunctions_lua.h"
#include "../PhysicsMaterial/MultiGroupXS/lua/multigroup_xs_lua_utils.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)
#define LUA_CMACRO1(x,y) \
        lua_pushnumber(L, y); \
        lua_setglobal(L, #x)

void chi_physics::lua_utils::RegisterLuaEntities(lua_State *L)
{
  //=================================== Solver
//  LUA_FMACRO1(chiSolverInitialize);
//  LUA_FMACRO1(chiSolverExecute);
//  LUA_FMACRO1(chiSolverStep);
//  LUA_FMACRO1(chiSolverSetBasicOption);
//  LUA_FMACRO1(chiSolverGetName);
//  LUA_FMACRO1(chiSolverGetFieldFunctionList);

//  lua_register(L, "chiSolverInitialize", temp::chiSolverInitialize);
//  lua_register(L, "chiSolverExecute", temp::chiSolverExecute);
//  lua_register(L, "chiSolverStep", temp::chiSolverStep);
//  lua_register(L, "chiSolverSetBasicOption", temp::chiSolverSetBasicOption);
//  lua_register(L, "chiSolverGetName", temp::chiSolverGetName);
//  lua_register(L, "chiSolverGetFieldFunctionList", temp::chiSolverGetFieldFunctionList);

  //=================================== Materials
  LUA_FMACRO1(chiPhysicsAddMaterial);
  LUA_FMACRO1(chiPhysicsMaterialAddProperty);
  LUA_FMACRO1(chiPhysicsMaterialSetProperty);
  LUA_FMACRO1(chiPhysicsMaterialGetProperty);
//  LUA_FMACRO1(chiPhysicsMaterialModifyTotalCrossSection);

  //=================================== Field functions
  chi_physics::field_function_lua_utils::RegisterLuaEntities(L);

  //=================================== Transport Cross-sections
  LUA_FMACRO1(chiPhysicsTransportXSCreate);
  LUA_FMACRO1(chiPhysicsTransportXSSet);
  LUA_FMACRO1(chiPhysicsTransportXSMakeCombined);
  LUA_FMACRO1(chiPhysicsTransportXSSetCombined);
  LUA_FMACRO1(chiPhysicsTransportXSGet);
  LUA_FMACRO1(chiPhysicsTransportXSExportToChiTechFormat);

  //Property indices
  LUA_CMACRO1(SCALAR_VALUE,           1);
  LUA_CMACRO1(TRANSPORT_XSECTIONS,    10);
  LUA_CMACRO1(ISOTROPIC_MG_SOURCE,    11);

  //Operation indices
  LUA_CMACRO1(SINGLE_VALUE,            0);
  LUA_CMACRO1(FROM_ARRAY,              1);
  LUA_CMACRO1(SIMPLEXS0,              20);
  LUA_CMACRO1(SIMPLEXS1,              21);
  LUA_CMACRO1(EXISTING,               22);
  LUA_CMACRO1(CHI_XSFILE,             23);

}