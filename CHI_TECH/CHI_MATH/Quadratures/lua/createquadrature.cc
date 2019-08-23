#include "../../../CHI_LUA/chi_lua.h"
#include "../../chi_math.h"
#include "../quadrature_gausslegendre.h"
#include "../quadrature_gausschebyshev.h"

extern CHI_MATH    chi_math_handler;


/**\defgroup LuaQuadrature Quadrature rules
 * \ingroup LuaMath*/

//########################################################## Create empty system
/** Creates a quadrature.
 *
\param QuadratureType int Quadrature identifier.
\param NumberOfPoints int Number of quadrature points.

Identifiers:\n
 GAUSS_LEGENDRE = Gauss-Legendre quadrature.
 GAUSS_CHEBYSHEV = Gauss-Chebyshev quadrature.

\return Returns a unique handle to the created quadrature rule

\ingroup LuaQuadrature
\author Jan*/
int chiCreateQuadrature(lua_State *L)
{

  //============================================= Parse argument
  int ident = lua_tonumber(L,1);
  int N = lua_tonumber(L,2);

  if (ident == 1) //GAUSS_LEGENDRE
  {
    printf("Creating Gauss-Legendre Quadrature\n");
    CHI_QUADRATURE_GAUSSLEGENDRE* new_quad = new CHI_QUADRATURE_GAUSSLEGENDRE;
    new_quad->Initialize(N,1000,1.0e-10,true);
    chi_math_handler.quadratures.push_back(new_quad);
    int index = chi_math_handler.quadratures.size()-1;
    lua_pushnumber(L,index);
    return 1;
  }
  else if (ident == 2) //GAUSS_CHEBYSHEV
  {
    printf("Creating Gauss-Chebyshev Quadrature\n");
    CHI_QUADRATURE_GAUSSCHEBYSHEV* new_quad = new CHI_QUADRATURE_GAUSSCHEBYSHEV;
    new_quad->Initialize(N,true);
    chi_math_handler.quadratures.push_back(new_quad);
    int index = chi_math_handler.quadratures.size()-1;
    lua_pushnumber(L,index);
    return 1;
  }
  return 0;
}