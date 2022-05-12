#include "ChiLua/chi_lua.h"

#include "chi_runtime.h"

#include "ChiMath/chi_math.h"
#include "../quadrature_gausslegendre.h"
#include "../quadrature_gausschebyshev.h"

#include "chi_log.h"
;

//########################################################## Create empty system
/** Creates a quadrature.
 *
\param QuadratureType int Quadrature identifier.
\param NumberOfPoints int Number of quadrature points.
\param VerboseFlag bool As the name implies. Default: false.

##_

###QuadratureType:\n
 GAUSS_LEGENDRE = Gauss-Legendre quadrature.
 GAUSS_CHEBYSHEV = Gauss-Chebyshev quadrature.

\return Returns a unique handle to the created quadrature rule

\ingroup LuaQuadrature
\author Jan*/
int chiCreateQuadrature(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (not ((num_args == 2) or (num_args == 3)))
    LuaPostArgAmountError("chiCreateQuadrature",2,num_args);

  LuaCheckNilValue("chiCreateQuadrature",L,1);
  LuaCheckNilValue("chiCreateQuadrature",L,2);

  //============================================= Parse argument
  int ident = lua_tonumber(L,1);
  int N = lua_tonumber(L,2);
  bool verbose = false;
  if (num_args == 3)
    verbose = lua_toboolean(L,3);

  if (ident == 1) //GAUSS_LEGENDRE
  {
    chi::log.Log() << "Creating Gauss-Legendre Quadrature\n";
    auto new_quad = std::make_shared<chi_math::QuadratureGaussLegendre>(N, verbose);
    chi::quadrature_stack.push_back(new_quad);
    int index = (int)chi::quadrature_stack.size()-1;
    lua_pushnumber(L,index);
    return 1;
  }
  else if (ident == 2) //GAUSS_CHEBYSHEV
  {
    chi::log.Log() << "Creating Gauss-Chebyshev Quadrature\n";
    auto new_quad = std::make_shared<chi_math::QuadratureGaussChebyshev>(N, verbose);
    chi::quadrature_stack.push_back(new_quad);
    int index = (int)chi::quadrature_stack.size()-1;
    lua_pushnumber(L,index);
    return 1;
  }
  return 0;
}