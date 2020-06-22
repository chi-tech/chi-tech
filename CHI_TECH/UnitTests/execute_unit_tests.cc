#include <ChiLua/chi_lua.h>

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//Prototypes
bool UnitTest_VectorNX();
bool UnitTest_MatrixNXxNX();

//###################################################################
/**Executes all the unit tests.
 */
int chiExecuteUnitTests(lua_State* L)
{
  chi_log.Log(LOG_0) << "Executing Chi-Tech unit tests.";

  //UnitTest_VectorNX();
  UnitTest_MatrixNXxNX();

  chi_log.Log(LOG_0) << "Chi-Tech unit testing completed.";
  return 0;
}