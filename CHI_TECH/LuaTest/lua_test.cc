#include <ChiLua/chi_lua.h>

#include "ChiMath/SparseMatrix/chi_math_sparse_matrix.h"

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**This is a lua test function.
\param argument1 Any Argument of any type.
\ingroup LuaGeneralUtilities
 */
int chiLuaTest(lua_State* L)
{
  chi_math::SparseMatrix A(5,5);

  std::cout << A.PrintS();

  A.Insert(0,4,20.0);
  A.Insert(0,0,2.0);
  A.Insert(0,1,-1.0);
  A.InsertAdd(0,1,-1.0);

  A.Compress();

  std::cout << std::endl;

  std::cout << A.PrintS();

  std::cout << std::endl;

  auto& indicesJ_rowI = A.rowI_indices[0];
  auto& valuesJ_rowI  = A.rowI_values[0];

  auto j = indicesJ_rowI.begin();
  auto v = valuesJ_rowI.begin();
  for (; j!= indicesJ_rowI.end(); ++j,++v)
  {
    std::cout << *j << " " << *v << "\n";
  }
  std::cout << std::endl;

  return 0;
}