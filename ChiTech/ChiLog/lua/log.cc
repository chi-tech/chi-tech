#include <chi_lua.h>

#include <chi_log.h>

#include "lua/chi_log_lua.h"

extern ChiLog& chi_log;

/** \defgroup LuaLogging D Output and Logging
 * \ingroup LuaUtilities*/

#define LUA_FMACRO1(x) lua_register(L, #x, x)
#define LUA_CMACRO1(x,y) \
        lua_pushnumber(L, y); \
        lua_setglobal(L, #x)

//###################################################################
/** Sets the verbosity level of the Logger.
 * This lua command will overwrite the currently set value.

\param int_level int Integer denoting verbosity level. Can be 0,1 or 2
 [default:0]

\ingroup LuaLogging
\author Jan*/
int chiLogSetVerbosity(lua_State *L)
{
  int num_args = lua_gettop(L);

  if (num_args == 0) {return 0;}
  else
  {
    int level = lua_tonumber(L,1);
    if (level<=2)
    {
      chi_log.SetVerbosity(level);
    }
  }
  return 0;
}

//###################################################################
/**Logs a message depending on the log type specified.

\param LogType int Can be any of the log types specified below.
\param LogMsg char Message or value to be output to the log.

##_

### LogType
LOG_0\n
Write a log only if location 0.\n

LOG_0WARNING\n
Write a log only if location 0 and format it as a warning.\n

LOG_0ERROR\n
Write a log only if location 0 and format it as an error.\n

LOG_0VERBOSE_0\n
Same as LOG_0.\n

LOG_0VERBOSE_1\n
Write a log only if location 0 and the verbosity level is greater
or equal to 1.\n

LOG_0VERBOSE_2\n
Write a log only if location 0 and the verbosity level is greater
or equal to 1.\n

LOG_ALL, LOG_ALLWARNING, LOG_ALLERROR,\n
LOG_ALLVERBOSE_0, LOG_ALLVERBOSE_1, LOG_ALLVERBOSE_2\n
Has the same meaning as their LOG_0 counterparts but instead applies to
all locations in the parallel context.

\ingroup LuaLogging
\author Jan
*/
int chiLog(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args != 2)
    LuaPostArgAmountError("chiLog",2,num_args);

  int         mode    = lua_tonumber(L,1);
  const char* message = lua_tostring(L,2);

  chi_log.Log((LOG_LVL)mode) << message << std::endl;

  return 0;
}

void chi_log_utils::lua_utils::RegisterLuaEntities(lua_State *L)
{
  LUA_FMACRO1(chiLogSetVerbosity);
  LUA_FMACRO1(chiLog);
  LUA_CMACRO1(LOG_0,          1);
  LUA_CMACRO1(LOG_0WARNING,   2);
  LUA_CMACRO1(LOG_0ERROR,     3);
  LUA_CMACRO1(LOG_0VERBOSE_0, 4);
  LUA_CMACRO1(LOG_0VERBOSE_1, 5);
  LUA_CMACRO1(LOG_0VERBOSE_2, 6);
  LUA_CMACRO1(LOG_ALL,          7);
  LUA_CMACRO1(LOG_ALLWARNING,   8);
  LUA_CMACRO1(LOG_ALLERROR,     9);
  LUA_CMACRO1(LOG_ALLVERBOSE_0, 10);
  LUA_CMACRO1(LOG_ALLVERBOSE_1, 11);
  LUA_CMACRO1(LOG_ALLVERBOSE_2, 12);
}