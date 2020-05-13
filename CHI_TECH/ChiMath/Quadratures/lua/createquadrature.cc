#include "../../../ChiLua/chi_lua.h"
#include "../../chi_math.h"
#include "../quadrature_gausslegendre.h"
#include "../quadrature_gausschebyshev.h"

extern ChiMath&     chi_math_handler;

#include <chi_log.h>

extern ChiLog& chi_log;


/**\defgroup LuaQuadrature Quadrature rules
 * \ingroup LuaMath*/

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
  size_t num_args = lua_gettop(L);

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
    chi_log.Log(LOG_0) << "Creating Gauss-Legendre Quadrature\n";
    chi_math::QuadratureGaussLegendre* new_quad = new chi_math::QuadratureGaussLegendre;
    new_quad->Initialize(N,1000,1.0e-10,verbose);
    chi_math_handler.quadratures.push_back(new_quad);
    int index = chi_math_handler.quadratures.size()-1;
    lua_pushnumber(L,index);
    return 1;
  }
  else if (ident == 2) //GAUSS_CHEBYSHEV
  {
    chi_log.Log(LOG_0) << "Creating Gauss-Chebyshev Quadrature\n";
    chi_math::QuadratureGaussChebyshev* new_quad = new chi_math::QuadratureGaussChebyshev;
    new_quad->Initialize(N,verbose);
    chi_math_handler.quadratures.push_back(new_quad);
    int index = chi_math_handler.quadratures.size()-1;
    lua_pushnumber(L,index);
    return 1;
  }
  return 0;
}