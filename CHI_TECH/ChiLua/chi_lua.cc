#include "chi_lua.h"

#include <string>
#include <sstream>

#include <chi_log.h>

extern ChiLog chi_log;

void LuaPostArgAmountError(const char* func_name,int expected, int given)
{
  chi_log.Log(LOG_ALLERROR)
  << "Incorrect amount of arguments supplied in "
  << std::string(func_name)
  << " expected " << expected << " arguments "
  << " but " << given << " provided";
  exit(EXIT_FAILURE);
}

void LuaCheckNilValue(const char* func_name, lua_State* L, int arg)
{
  if (lua_isnil(L,arg))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Nil -value supplied in "
      << std::string(func_name)
      << " argument " << arg;
    exit(EXIT_FAILURE);
  }
}

std::string LuaSourceInfo(lua_State* L, const char* func_name)
{
  lua_Debug err_info;
  lua_getstack(L,1,&err_info);
  lua_getinfo(L, "nSl", &err_info);

  std::stringstream ret_str;
  ret_str << func_name
  << " " << err_info.source
  << " line " << err_info.currentline;

  return ret_str.str();
}