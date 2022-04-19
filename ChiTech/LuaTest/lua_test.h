#ifndef CHITECH_LUA_TEST_H
#define CHITECH_LUA_TEST_H

#include "chi_lua.h"
int chiLuaTest(lua_State* L);

namespace chi_lua_test
{
  namespace lua_utils
  {
    void RegisterLuaEntities(lua_State* L);
  }//namespace lua_utils
}//namespace chi_lua_test

#endif //CHITECH_LUA_TEST_H
