#include "chi_lua.h"

#include "chi_runtime.h"

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
// NOLINTBEGIN(misc-no-recursion)
/**This function recursively processes table values. If the value is
 * a primitive type the recursion stops and the parameter block, which is
 * currently active, will be extended with a parameter of this primitive
 * type. If the value is another table, a new `Block`-type will be instantiated
 * and the table recursion will then operate on this new block.*/
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
      block.MakeAddParameter(key_str_name, bool_value);
      break;
    }
    case LUA_TNUMBER:
    {
      if (lua_isinteger(L, -1))
      {
        const int64_t number_value = lua_tointeger(L, -1);
        block.MakeAddParameter(key_str_name, number_value);
      }
      else
      {
        const double number_value = lua_tonumber(L, -1);
        block.MakeAddParameter(key_str_name, number_value);
      }

      break;
    }
    case LUA_TSTRING:
    {
      const std::string string_value = lua_tostring(L, -1);
      block.MakeAddParameter(key_str_name, string_value);
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
// NOLINTEND(misc-no-recursion)


//###################################################################
// NOLINTBEGIN(misc-no-recursion)
/**This function operates on table keys recursively. It has a specific
 * behavior if it detects an array.*/
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

      if (block.Type() != chi_data_types::ParameterBlockType::ARRAY)
        block.ChangeToArray();

      number_key_encountered = true;
      const std::string key_str_name = std::to_string(key_number_index);
      RecursivelyParseTableValues(L, block, key_str_name);
      ++key_number_index;
    }

    lua_pop(L, 1);
  }
}
// NOLINTEND(misc-no-recursion)

//###################################################################
/**This is the root command for parsing a table as a parameter block.
 * Example table:
\code
block =
{
  enabled = true,
  it_method = "gmres",
  nl_abs_tol = 1.0e-12,
  nl_max_its = 33,
  sub1 =
  {
    ax_method = 2,
    l_abs_tol = 1.0e-2
  },
  sub2 =
  {
    ax_method = 3,
    l_abs_tol = 1.0e-3,
    blocks = {99, 98, 97},
    cblocks = {{1,2,3},{4,5,6},{7,8,9}}
  }
}

chiUnitTests_Test_paramblock(--[[verbose=]]true, block)
\endcode*/
std::shared_ptr<chi_data_types::ParameterBlock> TableParserAsParameterBlock::
  ParseTable(lua_State* L, int table_stack_index)
{
  auto param_block = std::make_shared<ParamBlock>();

  RecursivelyParseTableKeys(L, table_stack_index, *param_block);

  return param_block;
}

}//namespace chi_lua