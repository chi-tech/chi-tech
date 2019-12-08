#include <ChiLua/chi_lua.h>

#include <chi_log.h>
extern ChiLog chi_log;

//###################################################################
/**This is a lua test function.
\param argument1 Any Argument of any type.
\ingroup LuaGeneralUtilities
 */
int chiLuaTest(lua_State* L)
{
  //============================== Optional argument protection
  int num_args = lua_gettop(L);
  if (num_args<1)
    LuaPostArgAmountError("chiLuaTest",1,num_args);

  //============================== Obtain first argument from stack
  const char* argument_1 = lua_tostring(L,1);

  //============================== Print to screen
  chi_log.Log(LOG_ALL) << "LuaTest: " << argument_1 << std::endl;

  size_t tag = chi_log.GetRepeatingEventTag(std::string());


  chi_log.LogEvent(tag,
                   ChiLog::EventType::SINGLE_OCCURRENCE,
                   std::make_shared<ChiLog::EventInfo>(std::string("A"),2.0));
  chi_log.LogEvent(tag,
                   ChiLog::EventType::SINGLE_OCCURRENCE,
                   std::make_shared<ChiLog::EventInfo>(std::string("B")));
  chi_log.LogEvent(tag,
                   ChiLog::EventType::SINGLE_OCCURRENCE,
                   std::make_shared<ChiLog::EventInfo>(std::string("C"),2.0));

  chi_log.Log(LOG_ALL)
    << chi_log.ProcessEvent(tag, ChiLog::EventOperation::AVERAGE_VALUE);
  std::cout << chi_log.PrintEventHistory(tag);


  return 0;
}