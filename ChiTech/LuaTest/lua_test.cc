#include <ChiLua/chi_lua.h>

#include "ChiMath/Quadratures/quadrature_triangle.h"

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
  int order_in = lua_tonumber(L,1);

  chi_log.Log() << "Hi";
  auto order = (chi_math::QuadratureOrder)order_in;

  chi_math::QuadratureTriangle qtri(order);

  {
    auto& quad = qtri;

    size_t np = quad.qpoints.size();
    for (size_t qp=0; qp<np; ++qp)
      chi_log.Log()
        << qp << " " << quad.weights[qp] << " " << quad.qpoints[qp].PrintS()
        << (quad.qpoints[qp].x + quad.qpoints[qp].y < 1.0);

    double weight_sum = 0.0;
    for (auto w : quad.weights) weight_sum += w;

    chi_log.Log() << weight_sum;
  }


  return 0;
}