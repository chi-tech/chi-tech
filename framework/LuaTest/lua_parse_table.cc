#include "chi_lua.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "ChiDataTypes/chi_data_types.h"
#include "ChiDataTypes/parameter_block.h"

#define MakeParamBlock std::make_unique<ParamBlock>

#define ExceptionLuaNilValue \
throw std::logic_error(std::string(__PRETTY_FUNCTION__) + \
": Encountered nil value assigned to key " + key_str_name)

#define ExceptionLuaUnsupportedValue \
throw std::logic_error(std::string(__PRETTY_FUNCTION__) + \
": Encountered unsupported value type " + lua_typename(L, lua_type(L, -2)) + \
" for key " + key_str_name)

#define ExceptionMixStringNumberKeys \
throw std::logic_error(std::string(__PRETTY_FUNCTION__) + \
": Encountered mixed key types (string and number)")

namespace chi_lua
{

typedef chi_data_types::ParameterBlock ParamBlock;


//###################################################################
/***/
void TableParserAsParameterBlock::
  RecursivelyParseTableValues(
    lua_State* L, ParamBlock &block, const std::string& key_str_name)
{
  switch (lua_type(L, -1))
  {
    case LUA_TNIL:     ExceptionLuaNilValue;
    case LUA_TBOOLEAN:
    {
      const bool bool_value = lua_toboolean(L, -1);
      block.AddParameter(MakeParamBlock(key_str_name, bool_value));
      break;
    }
    case LUA_TNUMBER:
    {
      if (lua_isinteger(L, -1))
      {
        const int64_t number_value = lua_tointeger(L, -1);
        block.AddParameter(MakeParamBlock(key_str_name, number_value));
      }
      else
      {
        const double number_value = lua_tonumber(L, -1);
        block.AddParameter(MakeParamBlock(key_str_name, number_value));
      }

      break;
    }
    case LUA_TSTRING:
    {
      const std::string string_value = lua_tostring(L, -1);
      chi::log.Log() << key_str_name + " -- " + string_value;
      block.AddParameter(MakeParamBlock(key_str_name, string_value));
      chi::log.Log() << block.GetParamValue<std::string>(key_str_name);
      break;
    }
    case LUA_TTABLE:
    {
      auto new_block = MakeParamBlock(key_str_name);
      RecursivelyParseTableKeys(L, lua_gettop(L), *new_block);
      block.AddParameter(std::move(new_block));
      break;
    }
    default:           ExceptionLuaUnsupportedValue;
  }//switch on value types
}



//###################################################################
void TableParserAsParameterBlock::
  RecursivelyParseTableKeys(
    lua_State* L, int t, chi_data_types::ParameterBlock& block)
{
  bool number_key_encountered = false;
  bool string_key_encountered = false;

  int key_number_index = 0;

  lua_pushnil(L); //first key
  while (lua_next(L, t) != 0) //pops the key, pushes next key and value
  {
    if (lua_type(L, -2) == LUA_TSTRING)
    {
      if (number_key_encountered) ExceptionMixStringNumberKeys;

      string_key_encountered = true;
      const std::string key_str_name = lua_tostring(L, -2);
      RecursivelyParseTableValues(L, block, key_str_name);
    }//if key is string

    // If the key is a number then the following apply:
    // - This must be an array of items
    // - All the keys in the table must be numbers
    if (lua_type(L, -2) == LUA_TNUMBER)
    {
      if (string_key_encountered) ExceptionMixStringNumberKeys;

      if (block.Type() != chi_data_types::ParameterBlockType::Array)
        block.ChangeToArray();

      number_key_encountered = true;
      const std::string key_str_name = std::to_string(key_number_index);
      RecursivelyParseTableValues(L, block, key_str_name);
      ++key_number_index;
    }

    lua_pop(L, 1);
  }
}

//###################################################################
std::shared_ptr<chi_data_types::ParameterBlock> TableParserAsParameterBlock::
  ParseTable(lua_State* L, int table_stack_index)
{
  auto param_block = std::make_shared<ParamBlock>();

  RecursivelyParseTableKeys(L, table_stack_index, *param_block);

  return param_block;
}

}//namespace chi_lua