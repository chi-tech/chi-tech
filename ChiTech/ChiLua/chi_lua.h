#ifndef CHI_LUA_H
#define CHI_LUA_H

extern "C"
{
#include<lua.h>
#include<lualib.h>
#include<lauxlib.h>
}

#include <typeinfo>
#include <string>
#include <vector>

void LuaPostArgAmountError(const std::string& func_name,int expected, int given);
void LuaCheckNilValue(const std::string& func_name, lua_State* L, int arg);
void LuaCheckStringValue(const std::string& func_name, lua_State* L, int arg);
void LuaCheckBoolValue(const std::string& func_name, lua_State* L, int arg);
void LuaCheckNumberValue(const std::string& func_name, lua_State* L, int arg);
void LuaCheckIntegerValue(const std::string& func_name, lua_State* L, int arg);
void LuaCheckTableValue(const std::string& func_name, lua_State* L, int arg);
std::string LuaSourceInfo(lua_State* L, const char* func_name);
void LuaPopulateVectorFrom1DArray(const std::string& func_name,
                                  lua_State* L,
                                  int table_arg_index,
                                  std::vector<double>& vec);

#endif