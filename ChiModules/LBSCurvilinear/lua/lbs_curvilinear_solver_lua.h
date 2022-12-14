#ifndef CHITECH_LBS_CURVILINEAR_SOLVER_LUA_H
#define CHITECH_LBS_CURVILINEAR_SOLVER_LUA_H

#include "chi_lua.h"

int chiLBSCurvilinearCreateSolver(lua_State *L);

namespace lbs_curvilinear::lua_utils
  {
    void RegisterLuaEntities(lua_State *L);
  }//namespace LinearBoltzmann

#endif //CHITECH_LBS_CURVILINEAR_SOLVER_LUA_H
