#include <ChiLua/chi_lua.h>
#include <iostream>

//###################################################################
///This is a lua test function.
///\param argument1 Any Argument of any type.
int chiLuaTest(lua_State* L)
{
  //============================== Optional argument protection
  int num_args = lua_gettop(L);
  if (num_args<1)
    LuaPostArgAmountError("chiLuaTest",1,num_args);

  //============================== Obtain first argument from stack
  const char* argument_1 = lua_tostring(L,1);

  //============================== Print to screen
  std::cout << "LuaTest: " << argument_1 << std::endl;

  return 0;
}