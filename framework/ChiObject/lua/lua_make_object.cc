#include "chi_lua.h"

#include "ChiObject/object_maker.h"
#include "ChiConsole/chi_console.h"

namespace chi_objects::lua_utils
{

int chiMakeObject(lua_State* L);

ChiConsoleRegisterLuaFunction(chiMakeObject);

// #############################################################################
/**Generic lua routine for the creation of objects.
 * \param type string. The type to create.
 * \param params ParameterBlock. A single block tree.

*/
int chiMakeObject(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckStringValue(fname, L, 1);
  LuaCheckTableValue(fname, L, 2);

  const std::string type = lua_tostring(L, 1);
  const auto params = chi_lua::TableParserAsParameterBlock::ParseTable(L, 2);

  const auto& object_maker = chi_objects::ObjectMaker::GetInstance();
  const size_t handle = object_maker.MakeObjectType(type, params);

  lua_pushinteger(L, static_cast<lua_Integer>(handle));
  return 1;
}

} // namespace chi_objects::lua_utils
