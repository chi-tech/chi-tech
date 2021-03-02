#include <ChiLua/chi_lua.h>

#include "ChiMath/dynamic_vector.h"
#include "ChiMath/dynamic_matrix.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**This is a lua test function.
\param argument1 Any Argument of any type.
\ingroup LuaGeneralUtilities
 */
int chiLuaTest(lua_State* L)
{

  chi_math::DynamicVector<double> vec(5, 1.0);
  chi_math::DynamicMatrix<double> mat(5,7,1.0);

  chi_log.Log() << vec.PrintStr() << std::endl;
  chi_log.Log() << "Hello\n" << mat.PrintStr();

  return 0;
}