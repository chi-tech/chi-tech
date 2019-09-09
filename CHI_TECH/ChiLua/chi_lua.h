#ifndef CHI_LUA_H
#define CHI_LUA_H

extern "C"
{
#include<lua.h>
#include<lualib.h>
#include<lauxlib.h>
}

#include <typeinfo>

void LuaPostArgAmountError(const char* func_name,int expected, int given);
void LuaCheckNilValue(const char* func_name, lua_State* L, int arg);

#endif